/*
	This program handles CODAR Time Series (TS) files, as ncdump and ncgen do for NetCDF data.
	The same executable can be called tsdump or tsgen to act as follows:
	- tsdump reads a binary TS file and generates an ASCII text representation of the data that can then be edited.
	- tsgen reads an ascii file produced by tsdump and converts it into a binary TS file.

	(c) 2018 Marcel Losekoot, Bodega Marine Laboratory, UC Davis.

	Notes: the binary TS file is bigendian by definition, so the program tests itself and corrects accordingly.
*/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>		// ?
#include <time.h>		// ctime()
#include <ctype.h>		// islower(), isprint()
#include <unistd.h>		// fseek constant SEEK_END
#include <string.h>		// memset()
#include <stdint.h>		// uint32_t
#include <libgen.h>		// basename()
#include <math.h>		// round()

typedef uint32_t fourcc ;	// four bytes that are subject to byte swapping

// declare a node for the linked list
struct node
{
	fourcc key ;		// the key name
	uint32_t size ;		// the size of the node's data block
	unsigned char *data ;	// the data block
	struct node *next ;	// the next node
} ;

// define key value codes, used to recognize and label block types
#define KEY_AQLV	(fourcc )0x4151564c	// "AQVL"
#define KEY_HEAD	(fourcc )0x48454144	// "HEAD"
#define KEY_sign	(fourcc )0x7369676e	// "sign"
#define KEY_mcda	(fourcc )0x6d636461	// "mcda"
#define KEY_cnst	(fourcc )0x636e7374	// "cnst"
#define KEY_swep	(fourcc )0x73776570	// "swep"
#define KEY_fbin	(fourcc )0x6662696e	// "fbin"
#define KEY_BODY	(fourcc )0x424f4459	// "BODY"
#define KEY_gtag	(fourcc )0x67746167	// "gtag"
#define KEY_atag	(fourcc )0x61746167	// "atag"
#define KEY_indx	(fourcc )0x696e6478	// "indx"
#define KEY_scal	(fourcc )0x7363616c	// "scal"
#define KEY_alvl	(fourcc )0x616c766c	// "alvl"
#define KEY_END 	(fourcc )0x454e4420	// "END "

// define size constants for character strings
#define SIZE_DESCRIPTION	64
#define SIZE_OWNERNAME		64
#define SIZE_COMMENT		64

// define codes for recognizing binary format and type options
#define BINFORMAT_CVIQ	(fourcc )0x63766971	// "cviq"
#define BINTYPE_FLT4	(fourcc )0x666c7434	// "flt4"
#define BINTYPE_FIX2	(fourcc )0x66697832	// "fix2"
#define BINTYPE_FIX3	(fourcc )0x66697833	// "fix3"
#define BINTYPE_FIX4	(fourcc )0x66697834	// "fix4"

// define a struct to hold all the useful stuff
struct config
{
	// file summary data
	unsigned int count_all ;
	unsigned int count_aqlv ;
	unsigned int count_head ;
	unsigned int count_sign ;
	unsigned int count_mcda ;
	unsigned int count_cnst ;
	unsigned int count_swep ;
	unsigned int count_fbin ;
	unsigned int count_body ;
	unsigned int count_gtag ;
	unsigned int count_atag ;
	unsigned int count_indx ;
	unsigned int count_scal ;
	unsigned int count_alvl ;
	unsigned int count_end ;
	unsigned int count_other ;
	unsigned int count_bad ;
	unsigned long count_samples ;
	// data for file signature
	fourcc version ;			// file version
	fourcc filetype ;			// file type
	fourcc sitecode; 			// site code
	unsigned int userflags ;		// user flags
	char description[SIZE_DESCRIPTION] ;	// description
	char ownername[SIZE_OWNERNAME] ;	// owner name
	char comment[SIZE_COMMENT] ;		// comment
	// data for file timestamp
	time_t timestamp ;			// unix time of file creation
	// data for size information
	int nchannels ;				// number of channels
	int nsweeps ;				// number of sweeps
	int nsamples ;				// number of samples per sweep
	int iqindicator ;			// iq indicator, 1=? 2=IQ
	// sweep information
	int samplespersweep ;			// number of samples per sweep/channels (normally 2048)
	double sweepstart ;			// sweep start frequency in Hertz
	double sweepbandwidth ;			// sweep bandwidth in Hertz
	double sweeprate ;			// sweep rate in Hertz
	int rangeoffset ;			// rangeoffset (not used)
	// sample binary format information
	fourcc bin_format ;			// the binary format: normally 'cviq'
	fourcc bin_type ;			// the type of binary data: one of 'flt4', 'fix2', 'fix3' or 'fix4', normally fix2
	// gtag value
	int gtag ;				// unknown
	// atag value
	int atag ;				// unknown
	// index value
	int index ; 				// sweep index
	int min_index ;				// sweep index min value
	int max_index ;				// sweep index max value
	// scaling information
	double scalar_one ;			// scaling value for I
	double scalar_two ;			// scaling value for Q
} ;

struct block_header
{
	fourcc key ;
	uint32_t size ;
} __attribute__((packed)) ;	// disable padding to make the header line up with the file data

struct block_functions						// this struct is used to relate a key name with a set of functions
{
	fourcc key ;						// a 4 byte block key
	int (*fixup)(struct node *) ;				// a pointer to a function that is called to perform endian fixup on the data block
	int (*make)(struct node *, struct config *, FILE *) ;	// a pointer to a function that is called to create a data block from text
	int (*dump)(struct node *, struct config *, FILE *) ;	// a pointer to a function that is called to produce text output from a data block
	int (*gen)(struct node *, FILE *) ;			// a pointer to a function that is called to write out a binary version of the block
} ;


int check_little_endian(void) ;
void usage_tsdump(char *) ;
void usage_tsgen(char *) ;
int tsdump(FILE *, FILE *, int) ;
int tsgen(FILE *, FILE *) ;
int read_binary_file(FILE *, unsigned long, unsigned char *) ;
int check_header(unsigned char *) ;
struct node *parse_file(unsigned char *, unsigned long) ;
int parse_block(struct node *, unsigned char *, unsigned long) ;
int superblock(fourcc) ;
void show_list(struct node *) ;
void free_all_nodes(struct node *) ;
void free_all_nodes_and_data(struct node *) ;
void hexdump(unsigned char *, int) ;
void endian_fixup(void *, int) ;
void swapcopy(unsigned char *, unsigned char *, int) ;
void swapcopy2(unsigned char *, unsigned char *) ;
void swapcopy4(unsigned char *, unsigned char *) ;
void swapcopy8(unsigned char *, unsigned char *) ;
char *strkey(fourcc) ;
int fixup_sizes(struct node *) ;
uint32_t calculate_body_size(struct node *) ;
uint32_t calculate_head_size(struct node *) ;
int set_block_size(struct node *, fourcc , uint32_t) ;
int dump_list(struct node *, FILE *, int) ;
void chomp(char *, int) ;
int read_parameter(FILE *, char [], void *) ;
int fixup_data(struct node *) ;
struct block_functions *find_block_functions(fourcc) ;
int ts_write(struct node *, FILE *) ;
int count_alvl_lines(FILE *) ;

// a set of functions that dump the contents of a specific type of block
int dump_block_aqlv(struct node *, struct config *, FILE *) ;
int dump_block_head(struct node *, struct config *, FILE *) ;
int dump_block_sign(struct node *, struct config *, FILE *) ;
int dump_block_mcda(struct node *, struct config *, FILE *) ;
int dump_block_cnst(struct node *, struct config *, FILE *) ;
int dump_block_swep(struct node *, struct config *, FILE *) ;
int dump_block_fbin(struct node *, struct config *, FILE *) ;
int dump_block_body(struct node *, struct config *, FILE *) ;
int dump_block_gtag(struct node *, struct config *, FILE *) ;
int dump_block_atag(struct node *, struct config *, FILE *) ;
int dump_block_indx(struct node *, struct config *, FILE *) ;
int dump_block_scal(struct node *, struct config *, FILE *) ;
int dump_block_alvl(struct node *, struct config *, FILE *) ;
int dump_block_end(struct node *, struct config *, FILE *) ;

// a set of functions that fixup the data block of a specific type of block
int fixup_data_aqlv(struct node *) ;
int fixup_data_head(struct node *) ;
int fixup_data_sign(struct node *) ;
int fixup_data_mcda(struct node *) ;
int fixup_data_cnst(struct node *) ;
int fixup_data_swep(struct node *) ;
int fixup_data_fbin(struct node *) ;
int fixup_data_body(struct node *) ;
int fixup_data_gtag(struct node *) ;
int fixup_data_atag(struct node *) ;
int fixup_data_indx(struct node *) ;
int fixup_data_scal(struct node *) ;
int fixup_data_alvl(struct node *) ;
int fixup_data_end(struct node *) ;

// a set of functions that create a node for a specific type of block
int make_node_aqlv(struct node *, struct config *config, FILE *) ;
int make_node_head(struct node *, struct config *config, FILE *) ;
int make_node_sign(struct node *, struct config *config, FILE *) ;
int make_node_mcda(struct node *, struct config *config, FILE *) ;
int make_node_cnst(struct node *, struct config *config, FILE *) ;
int make_node_swep(struct node *, struct config *config, FILE *) ;
int make_node_fbin(struct node *, struct config *config, FILE *) ;
int make_node_body(struct node *, struct config *config, FILE *) ;
int make_node_gtag(struct node *, struct config *config, FILE *) ;
int make_node_atag(struct node *, struct config *config, FILE *) ;
int make_node_indx(struct node *, struct config *config, FILE *) ;
int make_node_scal(struct node *, struct config *config, FILE *) ;
int make_node_alvl(struct node *, struct config *config, FILE *) ;
int make_node_end(struct node *, struct config *config, FILE *) ;

// a set of functions that generate binary file data for a specific type of block
int gen_block_aqlv(struct node *, FILE *) ;
int gen_block_head(struct node *, FILE *) ;
int gen_block_sign(struct node *, FILE *) ;
int gen_block_mcda(struct node *, FILE *) ;
int gen_block_cnst(struct node *, FILE *) ;
int gen_block_swep(struct node *, FILE *) ;
int gen_block_fbin(struct node *, FILE *) ;
int gen_block_body(struct node *, FILE *) ;
int gen_block_gtag(struct node *, FILE *) ;
int gen_block_atag(struct node *, FILE *) ;
int gen_block_indx(struct node *, FILE *) ;
int gen_block_scal(struct node *, FILE *) ;
int gen_block_alvl(struct node *, FILE *) ;
int gen_block_end(struct node *, FILE *) ;


int Global_flag_little_endian = 1 ;	// 1 indicates this code is little endian, 0 means it's big endian. The binary file is big endian.


int main(int argc, char *argv[])
{
	char *program_name = basename(argv[0]) ;
	Global_flag_little_endian = check_little_endian() ;
	int err = 0 ;
	FILE *fdin ;
	FILE *fdout ;
	if( strcmp(program_name,"tsdump") == 0 )
	{
		// do tsdump
		if( argc < 3 )
		{
			usage_tsdump(program_name) ;
			return 0 ;
		}
		int just_header = 0 ;
		if( strcmp(argv[1],"-h") == 0 )
		{
			just_header = 1 ;
			argv++ ;
			argc-- ;
		}
		char *infilename = argv[1] ;
		if( (fdin = fopen(infilename,"rb")) == NULL )
		{
			printf("Cannot open input file '%s'\n",infilename) ;
			return 1 ;
		}
		char *outfilename = argv[2] ;
		if( (fdout = fopen(outfilename,"wt")) == NULL )
		{
			printf("Cannot open output file '%s'\n",outfilename) ;
			fclose(fdin) ;
			return 1 ;
		}
		err = tsdump(fdin,fdout,just_header) ;
	}
	if( strcmp(program_name,"tsgen") == 0 )
	{
		// do tsgen
		if( argc < 3 )
		{
			usage_tsgen(program_name) ;
			return 0 ;
		}
		char *infilename = argv[1] ;
		if( (fdin = fopen(infilename,"rt")) == NULL )
		{
			printf("Cannot open input file '%s'\n",infilename) ;
			return 1 ;
		}
		char *outfilename = argv[2] ;
		if( (fdout = fopen(outfilename,"wb")) == NULL )
		{
			printf("Cannot open output file '%s'\n",outfilename) ;
			fclose(fdout) ;
			return 1 ;
		}
		err = tsgen(fdin,fdout) ;
	}
	fclose(fdin) ;
	fclose(fdout) ;
	return err ;
}

int check_little_endian(void)	// check whether this code is big or little endian, returns 1 for little.
{
	unsigned int i = 1 ;
	char *c = (char *)&i ;
	if( *c )		
		return 1 ;	// little endian
	return 0 ;
}

void usage_tsdump(char *name)
{
	printf("Usage: %s [-h] infile outfile\n",name) ;
	printf("Processes CODAR SeaSonde TimeSeries data files.\n") ;
	printf("Reads a binary infile and writes an ascii text version to outfile.\n") ;
}

void usage_tsgen(char *name)
{
	printf("Usage: %s infile outfile\n",name) ;
	printf("Processes CODAR SeaSonde TimeSeries data files.\n") ;
	printf("Reads an ascii text infile and writes a binary version to outfile.\n") ;
}

int tsdump(FILE *infile, FILE *outfile, int just_header)
{
	fseek(infile,0L,SEEK_END) ;	// seek to the end of the file
	unsigned long filesize = ftell(infile) ;
	rewind(infile) ;
	unsigned char *filedata = malloc(filesize) ;	// try to make a buffer sized to read the entire file
	if( filedata == NULL )
	{
		printf("Cannot get memory for file with %ld bytes\n",filesize) ;
		return 1 ;
	}
	int err = 0 ;
	if( read_binary_file(infile,filesize,filedata) == 0 )
	{
		struct node *list = parse_file(filedata,filesize) ;
		if( list != NULL )
		{
			err = dump_list(list,outfile,just_header) ;
			free_all_nodes(list) ;
		}
	}
	free(filedata) ;
	return err ;
}

#define SIZE_LINE 80

int tsgen(FILE *infile, FILE *outfile)
{
	char line[SIZE_LINE] ;
	long line_count = 0 ;
	struct config config ;
	memset(&config,0,sizeof(struct config)) ;
	struct node root ;
	struct node *list = &root ;
	memset(list,0,sizeof(struct node)) ;
	while( fgets(line,SIZE_LINE,infile) )
	{
		chomp(line,SIZE_LINE) ;				// remove newline
		//printf("debug: chomped line is '%s'\n",line) ;
		line_count++ ;
		if( strlen(line) <= 1 ) continue ;		// skip empty lines
		if( index(line,':') != NULL ) continue ;	// skip parameter lines
		fourcc key = *(fourcc *)line ;			// extract the block type
		endian_fixup(&key,sizeof(key)) ;
		struct block_functions *block_functions = find_block_functions(key) ;	// returns a set of functions from the Global_functions_dictionary for this block type
		if( block_functions == NULL )
		{
			printf("Cannot gen block '%s'\n",strkey(key)) ;
			return 1 ;
		}
		int (*make_function)(struct node *, struct config *, FILE *) = block_functions->make ;
		int err = (*make_function)(list,&config,infile) ;	// calls the 'make' function from Function_dictionary corresponding to the block type
		if( err )
		{
			printf("Error in '%s' block starting at line %ld\n",strkey(key),line_count) ;
			return 1 ;
		}
		if( list->next != NULL ) list = list->next ;	// advance the list pointer to the newly created node
	}
	printf("Read %ld lines\n",line_count) ;
	fixup_sizes(&root) ;	// calculate body, head and aqlv block sizes, update nodes
	// write to outfile
	int err = ts_write(root.next,outfile) ;
	free_all_nodes_and_data(list) ;
	return err ;
}

int ts_write(struct node *list, FILE *outfile)
{
	//printf("debug: ts_write: start\n") ;
	while( list != NULL )
	{
		fourcc key = list->key ;
		struct block_functions *block_functions = find_block_functions(key) ;	// gets a set of functions from Global_functions_dictionary for this block type
		if( block_functions == NULL )
		{
			printf("Cannot write block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		int (*gen_function)(struct node *, FILE *) = block_functions->gen ;
		int err = (*gen_function)(list,outfile) ;	// calls the 'gen' function corresponding to the block type
		if( err )
		{
			printf("Error in '%s' block\n",strkey(key)) ;
			return 1 ;
		}
		list = list->next ;
	}
	//printf("debug: ts_write: finish\n") ;
	return 0 ;
}

void chomp(char *line, int max)		// delete trailing newline character
{
	line[max-1] = '\0' ;
	char *nl = index(line,0x0a) ;
	if( nl ) *nl = '\0' ;
}

int read_parameter(FILE *fd, char format[], void *buffer)	// looks for a line with format, copies the value to buffer
{
	size_t block_start = ftell(fd) ;	// put a finger in the file at the first line of the block
	int err = 0 ;
	char line[SIZE_LINE] ;
	//printf("debug: read_parameters: format='%s'\n",format) ;
	while( fgets(line,SIZE_LINE,fd) )
	{
		chomp(line,SIZE_LINE) ;
		if( strlen(line) == 0 )
		{
			printf("Cannot find parameter '%s'\n",format) ;
			err = 1 ;
			break ;			// stop on blank line
		}
		//printf("debug: read_parameters: line='%s'\n",line) ;
		int count = sscanf(line,format,buffer) ;
		//printf("debug: read_parameters: sscanf returned %d\n",count) ;
		if( count != 0 )
		{
			//printf("debug: read_parameters: using line %s\n",line) ;
			err = 0 ;
			break ;
		}
	}
	fseek(fd,block_start,SEEK_SET) ;
	return err ;
}

int read_binary_file(FILE * tsfile, unsigned long filesize, unsigned char *buffer)
{
	unsigned long count = fread(buffer,1,filesize,tsfile) ;
	if( count != filesize )
	{
		printf("Error reading ts file, only read %lu bytes out of %lu\n",count,filesize) ;
		return 1 ;
	}
	if( count > sizeof(struct block_header) )
	{
		if( check_header(buffer) )
			return 1 ;
	}
	return 0 ;
}

int check_header(unsigned char *buffer)
{
	struct block_header *header = (struct block_header *)buffer ;
	fourcc key = header->key ;
	endian_fixup(&key,sizeof(key)) ;
	//printf("debug: check_header: read %x as %x\n",header->key,key) ;
	if( key != KEY_AQLV )
	{
		printf("Bad header key: %x\n",key) ;
		return 1 ;
	}
	// TODO check file size (meh?)
	return 0 ;
}

struct node *parse_file(unsigned char *buffer, unsigned long length)
{
	struct node dummy_root ;	// use a dummy root node to prime the parser
	if( parse_block(&dummy_root,buffer,length) )
	{
		printf("Parser error\n") ;
		return NULL ;
	}
	return dummy_root.next ;	// return the next node as the true root
}

int parse_block(struct node *root, unsigned char *buffer, unsigned long length)
{
	while( length > 0 )
	{
		// make a new node to describe this block
		struct node *newnode = malloc(sizeof(struct node)) ;
		if( newnode == NULL )
		{
			printf("Malloc error\n") ;
			return 1 ;
		}
		memset(newnode,0,sizeof(struct node)) ;
		root->next = newnode ;					// link the previous node to the new node
		struct block_header *header = (struct block_header *)buffer ;	// make buffer contents accessible through struct block_header
		endian_fixup(&(header->key),sizeof(newnode->key)) ;	// fixup the endian order
		newnode->key = header->key ;				// copy the block type, aka key, from the header
		endian_fixup(&(header->size),sizeof(newnode->size)) ;	// fixup the endian order
		newnode->size = header->size ;				// copy the size of the data block (excluding header)
		length -= sizeof(struct block_header) ;			// reduce the block length by the size of the header
		buffer += sizeof(struct block_header) ;			// advance the buffer pointer by the size of the header
		if( newnode->size > length )
		{
			printf("Block '%s' size truncted from %u to %lu bytes\n",strkey(newnode->key),newnode->size,length) ;
			newnode->size = length ;
		}
		newnode->data = buffer ;				// point at the data portion of the block
		if( superblock(newnode->key) )				// if the block is a superblock, recursively parse its data block
		{
			if( parse_block(newnode,newnode->data,newnode->size) )
				return 1 ;
		}
		else
		{
			if( fixup_data(newnode) )			// otherwise, do endian fixup on the block's data, other stuff too (?)
				return 1 ;
		}
		// move on to the next block in the buffer
		length -= newnode->size ;				// reduce the block length by the size of the data block
		buffer += newnode->size ;				// advance the buffer pointer by the size of the data block
		do{							// iterate to deal with recursive parsing
			root = root->next ;				// advance the root pointer to the new node, ready to iterate
		}while( root->next != NULL ) ;
	}
	return 0 ;
}

struct block_functions Global_function_dictionary[] =		// a dictionary of  names and associated functions, used to lookup which function to call
{
	{ KEY_alvl, fixup_data_alvl, make_node_alvl, dump_block_alvl, gen_block_alvl  },
	{ KEY_gtag, fixup_data_gtag, make_node_gtag, dump_block_gtag, gen_block_gtag  },
	{ KEY_atag, fixup_data_atag, make_node_atag, dump_block_atag, gen_block_atag  },
	{ KEY_indx, fixup_data_indx, make_node_indx, dump_block_indx, gen_block_indx  },
	{ KEY_scal, fixup_data_scal, make_node_scal, dump_block_scal, gen_block_scal  },
	{ KEY_sign, fixup_data_sign, make_node_sign, dump_block_sign, gen_block_sign  },
	{ KEY_mcda, fixup_data_mcda, make_node_mcda, dump_block_mcda, gen_block_mcda  },
	{ KEY_cnst, fixup_data_cnst, make_node_cnst, dump_block_cnst, gen_block_cnst  },
	{ KEY_swep, fixup_data_swep, make_node_swep, dump_block_swep, gen_block_swep  },
	{ KEY_fbin, fixup_data_fbin, make_node_fbin, dump_block_fbin, gen_block_fbin  },
	{ KEY_AQLV, fixup_data_aqlv, make_node_aqlv, dump_block_aqlv, gen_block_aqlv  },
	{ KEY_HEAD, fixup_data_head, make_node_head, dump_block_head, gen_block_head  },
	{ KEY_BODY, fixup_data_body, make_node_body, dump_block_body, gen_block_body  },
	{ KEY_END, fixup_data_end, make_node_end, dump_block_end, gen_block_end  },
	{ 0, NULL, NULL, NULL, NULL }
} ;

struct block_functions *find_block_functions(fourcc key)
{
	if( key == 0 )
	{
		printf("Bad key (zero!)\n") ;
		return NULL ;
	}
	for( struct block_functions *block_functions = Global_function_dictionary ; block_functions->key != 0 ; block_functions++ )
	{
		if( key == block_functions->key )
			return block_functions ;
	}
	printf("Cannot locate functions for key '%s'\n",strkey(key)) ;
	return NULL ;
}

int fixup_data(struct node *node)	// figures out what type of block and what to do with it
{
	struct block_functions *block_functions = find_block_functions(node->key) ;	// returns a set of functions from the Global_functions_dictionary for this block type
	if( block_functions == NULL )
	{
		printf("Cannot fixup block '%s'\n",strkey(node->key)) ;
		return 1 ;
	}
	int (*fixup_function)(struct node *) = block_functions->fixup ;
	int err = (*fixup_function)(node) ;	// calls the fixup function corresponding to the block key
	if( err )
	{
		printf("Error fixing block %s\n",strkey(node->key)) ;
		return 1 ;
	}
	return 0 ;
}

void endian_fixup(void *original, int size)
{
	if( Global_flag_little_endian )
	{
		unsigned char tmp[8] ;
		memcpy(tmp,original,size) ;
		swapcopy(original,tmp,size) ;
	}
}

void swapcopy(unsigned char *dest, unsigned char *source, int size)
{
	switch( size )
	{
		case 2:
			swapcopy2(dest,source) ;
		break ;
		case 4:
			swapcopy4(dest,source) ;
		break ;
		case 8:
			swapcopy8(dest,source) ;
		break ;
	}
}

void swapcopy2(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[1] ;
	dest[1] = source[0] ;
}

void swapcopy4(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[3] ;
	dest[1] = source[2] ;
	dest[2] = source[1] ;
	dest[3] = source[0] ;
}

void swapcopy8(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[7] ;
	dest[1] = source[6] ;
	dest[2] = source[5] ;
	dest[3] = source[4] ;
	dest[4] = source[3] ;
	dest[5] = source[2] ;
	dest[6] = source[1] ;
	dest[7] = source[0] ;
}

int superblock(fourcc key)	// returns 1 if the key denotes a 'superblock', i.e. one that is composed of sub-blocks
{
	// check for one of the known superkeys
	if( key == KEY_AQLV ) return 1 ;
	if( key == KEY_HEAD ) return 1 ;
	if( key == KEY_BODY ) return 1 ;
	if( key == KEY_END ) return 1 ;
	return 0 ;		// no, this is not a superblock key
}

void hexdump(unsigned char *buffer, int n)
{
	for( int row = 0 ; row <= n/8 ; row++ )
	{
		for( int col = 0 ; col < 8 ; col++ )
		{
			int offset = row*8+col ;
			if( offset < n )
			{
				printf("%02x ",buffer[offset]) ;
			}
			else
			{
				printf("   ") ;
			}
		}
		printf("\t") ;
		for( int col = 0 ; col < 8 ; col++ )
		{
			int offset = row*8+col ;
			if( offset < n )
			{
				if( isprint(buffer[offset]) )
					printf("%c",buffer[offset]) ;
				else
					printf(".") ;
			}
			else
			{
				printf("   ") ;
			}
		}
		printf("\n") ;
	}
}

#define MAX_NAME 10
char Global_name[MAX_NAME] ;
char *strkey(fourcc name)
{
	int size = sizeof(name) ;
	if( size >= MAX_NAME )
		size = MAX_NAME-1 ;
	memcpy(Global_name,(void *)&name,size) ;
	endian_fixup(Global_name,size) ;
	Global_name[size] = '\0' ;
	return Global_name ;
}

int fixup_sizes(struct node *list)
{
	uint32_t head_size = calculate_head_size(list) ;
	if( head_size == 0 ) return 1 ;
	uint32_t body_size = calculate_body_size(list) ;
	if( body_size == 0 ) return 1 ;
	uint32_t aqlv_size = head_size + sizeof(struct block_header) + body_size + sizeof(struct block_header) ;
	if( set_block_size(list,KEY_AQLV,aqlv_size) ) return 1 ;
	if( set_block_size(list,KEY_HEAD,head_size) ) return 1 ;
	if( set_block_size(list,KEY_BODY,body_size) ) return 1 ;
	return 0 ;
}

uint32_t calculate_body_size(struct node *list)			// returns the size of the blocks in the BODY superblock
{
	int in_body = 0 ;
	uint32_t size = 0 ;
	while( list != NULL )
	{
		if( list->key == KEY_END )
			in_body = 0 ;
		if( in_body )
			size += ( list->size + sizeof(struct block_header) ) ;	// remember to count the block header
		if( list->key == KEY_BODY )
			in_body = 1 ;
		list = list->next ;
	} ;
	return size ;
}

uint32_t calculate_head_size(struct node *list)			// returns the size of the blocks in the HEAD superblock
{
	int in_head = 0 ;
	uint32_t size = 0 ;
	while( list != NULL )
	{
		if( list->key == KEY_END )
			in_head = 0 ;
		if( list->key == KEY_BODY )
			in_head = 0 ;
		if( in_head )
			size += (list->size + sizeof(struct block_header) ) ;	// remember to count the block header
		if( list->key == KEY_HEAD )
			in_head = 1 ;
		list = list->next ;
	} ;
	return size ;
}

int set_block_size(struct node *list, fourcc key, uint32_t size)	// returns 1 for failure, 0 for success
{
	while( list != NULL )
	{
		if( list->key == key )
		{
			list->size = size ;
			return 0 ;
		}
		list = list->next ;
	}
	return 1 ;
}

void show_list(struct node *list)
{
	unsigned int count = 0 ;
	while( list != NULL )
	{
		printf("Node %u: key %.4s size %u\n",count,(char *)&(list->key),list->size) ;
		list = list->next ;
		count++ ;
	}
}

int dump_list(struct node *list, FILE *outfile, int just_header) // goes through the list of nodes, writing an ascii text description of each node to outfile
{
	unsigned int count = 0 ;
	struct config config ;
	memset(&config,0,sizeof(struct config)) ;
	while( list != NULL )
	{
		//printf("debug: dump_list: node has key '%s'\n",strkey(list->key)) ;
		struct block_functions *block_functions = find_block_functions(list->key) ;	// returns a set of functions from the Global_functions_dictionary for this block type
		if( block_functions == NULL )
		{
			printf("Cannot dump block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		int (*dump_function)(struct node *, struct config *, FILE *) = block_functions->dump ;
		if( just_header && dump_function == dump_block_body ) return 0 ;
		int err = (*dump_function)(list,&config,outfile) ;	// calls a function from the Global_function_dictionary corresponding to the block key
		if( err )
		{
			printf("Error dumping block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		list = list->next ;
	}
	return 0 ;
}

int fixup_data_aqlv(struct node *node)
{
	return 0 ;
}

int dump_block_aqlv(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_AQLV)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_aqlv(struct node *list, struct config *config, FILE *fd)		// creates a new node for AQLV block, fd not used
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_AQLV ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_aqlv(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}

int fixup_data_head(struct node *node)
{
	return 0 ;
}

int dump_block_head(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_HEAD)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_head(struct node *list, struct config *config, FILE *fd)
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_HEAD ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_head(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_sign				// contains only the data portion for this type of block, i.e. no header data (key, size)
{
	fourcc version ;			// file version
	fourcc filetype ;			// file type
	fourcc sitecode ; 			// site code
	uint32_t userflags ;			// user flags
	char description[SIZE_DESCRIPTION] ;	// file description
	char ownername[SIZE_OWNERNAME] ;	// owner name
	char comment[SIZE_COMMENT] ;		// comment
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_sign(struct node *node)
{
	if( node->size < sizeof(struct block_sign) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_sign)) ;
		return 1 ;
	}
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	endian_fixup(&(sign->version),sizeof(sign->version)) ;
	endian_fixup(&(sign->filetype),sizeof(sign->filetype)) ;
	endian_fixup(&(sign->sitecode),sizeof(sign->sitecode)) ;
	endian_fixup(&(sign->userflags),sizeof(sign->userflags)) ;
	return 0 ;
}

int dump_block_sign(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_sign) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_sign)) ;
		return 1 ;
	}
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_sign)) ;
	fprintf(outfile,"version:%s\n",strkey(sign->version)) ;
	fprintf(outfile,"filetype:%s\n",strkey(sign->filetype)) ;
	fprintf(outfile,"sitecode:%s\n",strkey(sign->sitecode)) ;
	fprintf(outfile,"userflags:%x\n",sign->userflags) ;
	fprintf(outfile,"description:%.*s\n",SIZE_DESCRIPTION,sign->description) ;
	fprintf(outfile,"ownername:%.*s\n",SIZE_OWNERNAME,sign->ownername) ;
	fprintf(outfile,"comment:%.*s\n",SIZE_OWNERNAME,sign->comment) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_sign(struct node *list, struct config *config, FILE *fd)
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_sign ;
	struct block_sign *sign = malloc(sizeof(struct block_sign)) ;
	if( sign == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(sign,0,sizeof(struct block_sign)) ;
	newnode->data = (unsigned char *)sign ;
	newnode->size = sizeof(struct block_sign) ;
	if( read_parameter(fd,"version:%4c",(void *)&(sign->version)) ) return 1 ;
	if( read_parameter(fd,"filetype:%4c",(void *)&(sign->filetype)) ) return 1 ;
	if( read_parameter(fd,"sitecode:%4c",(void *)&(sign->sitecode)) ) return 1 ;
	if( read_parameter(fd,"userflags:%x",(void *)&(sign->userflags)) ) return 1 ;
	char format[32] ;
	sprintf(format,"description:%%%dc",SIZE_DESCRIPTION) ;
	if( read_parameter(fd,format,(void *)(sign->description)) ) return 1 ;
	sprintf(format,"ownername:%%%dc",SIZE_OWNERNAME) ;
	if( read_parameter(fd,format,(void *)(sign->ownername)) ) return 1 ;
	sprintf(format,"comment:%%%dc",SIZE_COMMENT) ;
	if( read_parameter(fd,format,(void *)(sign->comment)) ) return 1 ;
	return 0 ;
}

int gen_block_sign(struct node *node, FILE *outfile)
{
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->version),sizeof(sign->version)) ;
	if( fwrite(&(sign->version),sizeof(sign->version),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->filetype),sizeof(sign->filetype)) ;
	if( fwrite(&(sign->filetype),sizeof(sign->filetype),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->sitecode),sizeof(sign->sitecode)) ;
	if( fwrite(&(sign->sitecode),sizeof(sign->sitecode),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->userflags),sizeof(sign->userflags)) ;
	if( fwrite(&(sign->userflags),sizeof(sign->userflags),1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->description),SIZE_DESCRIPTION,1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->ownername),SIZE_OWNERNAME,1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->comment),SIZE_COMMENT,1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_mcda
{
	uint32_t timestamp ;		// mac timestamp of first sweep
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_mcda(struct node *node)
{
	if( node->size < sizeof(struct block_mcda) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_mcda)) ;
		return 1 ;
	}
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	endian_fixup(&(mcda->timestamp),sizeof(mcda->timestamp)) ;
	return 0 ;
}

int dump_block_mcda(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_mcda) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_mcda)) ;
		return 1 ;
	}
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	uint32_t mcda_time = mcda->timestamp ;
	time_t t = mcda_time ;
	fprintf(outfile,"%s\n",strkey(KEY_mcda)) ;
	if( t != 0 )
	{
		t -= 2082844800 ;	// move epoc from 1904-01-01 00:00:00 to 1970-01-01 00:00:00 
		fprintf(outfile,"timestamp:%lu (NB: seconds since 1970) (%.24s)\n",t,ctime(&t)) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_mcda(struct node *list, struct config *config, FILE *fd)		// creates a new node for mcda block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_mcda ;
	struct block_mcda *mcda = malloc(sizeof(struct block_mcda)) ;
	if( mcda == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(mcda,0,sizeof(struct block_mcda)) ;
	newnode->data = (unsigned char *)mcda ;
	newnode->size = sizeof(struct block_mcda) ;
	if( read_parameter(fd,"timestamp:%u",(void *)&(mcda->timestamp)) ) return 1 ;
	mcda->timestamp += 2082844800 ;	// move epoc from 1970-01-01 00:00:00 to 1904-01-01 00:00:00 
	return 0 ;
}

int gen_block_mcda(struct node *node, FILE *outfile)
{
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(mcda->timestamp),sizeof(mcda->timestamp)) ;
	if( fwrite(&(mcda->timestamp),sizeof(mcda->timestamp),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_cnst
{
	int32_t nchannels ;		// number of antennas/channels (normally 3)
	int32_t nsweeps ;		// number of sweeps (normally 32)
	int32_t nsamples ;		// number of samples (normally 2048)
	int32_t iqindicator ;		// iqindicator: 1=?, 2=IQ
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_cnst(struct node *node)
{
	if( node->size < sizeof(struct block_cnst) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_cnst)) ;
		return 1 ;
	}
	struct block_cnst *cnst = (struct block_cnst *)node->data ;
	endian_fixup(&(cnst->nchannels),sizeof(cnst->nchannels)) ;
	endian_fixup(&(cnst->nsweeps),sizeof(cnst->nchannels)) ;
	endian_fixup(&(cnst->nsamples),sizeof(cnst->nsamples)) ;
	endian_fixup(&(cnst->iqindicator),sizeof(cnst->iqindicator)) ;
	return 0 ;
}

int dump_block_cnst(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_cnst) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_cnst)) ;
		return 1 ;
	}
	struct block_cnst *cnst = (struct block_cnst *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_cnst)) ;
	fprintf(outfile,"nchannels:%d\n",cnst->nchannels) ;
	fprintf(outfile,"nsweeps:%d\n",cnst->nsweeps) ;
	fprintf(outfile,"nsamples:%d\n",cnst->nsamples) ;
	fprintf(outfile,"iqindicator:%d\n",cnst->iqindicator) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_cnst(struct node *list, struct config *config, FILE *fd)		// creates a new node for cnst block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_cnst ;
	struct block_cnst *cnst = malloc(sizeof(struct block_cnst)) ;
	if( cnst == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(cnst,0,sizeof(struct block_cnst)) ;
	newnode->data = (unsigned char *)cnst ;
	newnode->size = sizeof(struct block_cnst) ;
	if( read_parameter(fd,"nchannels:%d",(void *)&(cnst->nchannels)) ) return 1 ;
	if( read_parameter(fd,"nsweeps:%d",(void *)&(cnst->nsweeps)) ) return 1 ;
	if( read_parameter(fd,"nsamples:%d",(void *)&(cnst->nsamples)) ) return 1 ;
	if( read_parameter(fd,"iqindicator:%d",(void *)&(cnst->iqindicator)) ) return 1 ;
	return 0 ;
}

int gen_block_cnst(struct node *node, FILE *outfile)
{
	struct block_cnst *cnst = (struct block_cnst *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nchannels),sizeof(cnst->nchannels)) ;
	if( fwrite(&(cnst->nchannels),sizeof(cnst->nchannels),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nsweeps),sizeof(cnst->nsweeps)) ;
	if( fwrite(&(cnst->nsweeps),sizeof(cnst->nsweeps),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nsamples),sizeof(cnst->nsamples)) ;
	if( fwrite(&(cnst->nsamples),sizeof(cnst->nsamples),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->iqindicator),sizeof(cnst->iqindicator)) ;
	if( fwrite(&(cnst->iqindicator),sizeof(cnst->iqindicator),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_swep
{
	int32_t samplespersweep ;	// number of samples per sweep/channels (normally 2048)
	double sweepstart ;		// sweep start frequency in Hertz
	double sweepbandwidth ;		// sweep bandwidth in Hertz
	double sweeprate ;		// sweep rate in Hertz
	int32_t rangeoffset ;		// rangeoffset (not used)
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_swep(struct node *node)
{
	if( node->size < sizeof(struct block_swep) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_swep)) ;
		return 1 ;
	}
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	endian_fixup(&(swep->samplespersweep),sizeof(swep->samplespersweep)) ;
	endian_fixup(&(swep->sweepstart),sizeof(swep->sweepstart)) ;
	endian_fixup(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth)) ;
	endian_fixup(&(swep->sweeprate),sizeof(swep->sweeprate)) ;
	endian_fixup(&(swep->rangeoffset),sizeof(swep->rangeoffset)) ;
	return 0 ;
}

int dump_block_swep(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_swep) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_swep)) ;
		return 1 ;
	}
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_swep)) ;
	fprintf(outfile,"samplespersweep:%d\n",swep->samplespersweep) ;
	fprintf(outfile,"sweepstart:%.20lf\n",swep->sweepstart) ;
	fprintf(outfile,"sweepbandwidth:%.20lf\n",swep->sweepbandwidth) ;
	fprintf(outfile,"sweeprate:%.20lf\n",swep->sweeprate) ;
	fprintf(outfile,"rangeoffset:%d\n",swep->rangeoffset) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_swep(struct node *list, struct config *config, FILE *fd)		// creates a new node for swep block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_swep ;
	struct block_swep *swep = malloc(sizeof(struct block_swep)) ;
	if( swep == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(swep,0,sizeof(struct block_swep)) ;
	newnode->data = (unsigned char *)swep ;
	newnode->size = sizeof(struct block_swep) ;
	if( read_parameter(fd,"samplespersweep:%d",(void *)&(swep->samplespersweep)) ) return 1 ;
	if( read_parameter(fd,"sweepstart:%lf",(void *)&(swep->sweepstart)) ) return 1 ;
	if( read_parameter(fd,"sweepbandwidth:%lf",(void *)&(swep->sweepbandwidth)) ) return 1 ;
	if( read_parameter(fd,"sweeprate:%lf",(void *)&(swep->sweeprate)) ) return 1 ;
	if( read_parameter(fd,"rangeoffset:%d",(void *)&(swep->rangeoffset)) ) return 1 ;
	return 0 ;
}

int gen_block_swep(struct node *node, FILE *outfile)
{
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->samplespersweep),sizeof(swep->samplespersweep)) ;
	if( fwrite(&(swep->samplespersweep),sizeof(swep->samplespersweep),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweepstart),sizeof(swep->sweepstart)) ;
	if( fwrite(&(swep->sweepstart),sizeof(swep->sweepstart),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth)) ;
	if( fwrite(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweeprate),sizeof(swep->sweeprate)) ;
	if( fwrite(&(swep->sweeprate),sizeof(swep->sweeprate),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->rangeoffset),sizeof(swep->rangeoffset)) ;
	if( fwrite(&(swep->rangeoffset),sizeof(swep->rangeoffset),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_fbin
{
	fourcc bin_format ;	// format (normally 'cviq')
	fourcc bin_type ;	// type of ALVL data ('flt4','fix4','fix3','fix2')
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_fbin(struct node *node)
{
	if( node->size < sizeof(struct block_fbin) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_fbin)) ;
		return 1 ;
	}
	struct block_fbin *fbin = (struct block_fbin *)node->data ;
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;
	return 0 ;
}

int dump_block_fbin(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_fbin) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_fbin)) ;
		return 1 ;
	}
	struct block_fbin *fbin = (struct block_fbin *)(node->data) ;
	config->bin_type = fbin->bin_type ;	// remember this for alvl blocks
	fprintf(outfile,"%s\n",strkey(KEY_fbin)) ;
	fprintf(outfile,"format:%s\n",strkey(fbin->bin_format)) ;
	fprintf(outfile,"type:%s\n",strkey(fbin->bin_type)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_fbin(struct node *list, struct config *config, FILE *fd)		// creates a new node for fbin block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_fbin ;
	struct block_fbin *fbin = malloc(sizeof(struct block_fbin)) ;
	if( fbin == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(fbin,0,sizeof(struct block_fbin)) ;
	newnode->data = (unsigned char *)fbin ;
	newnode->size = sizeof(struct block_fbin) ;
	char format[16] ;
	sprintf(format,"format:%%%lus",sizeof(fbin->bin_format)) ;
	if( read_parameter(fd,format,(void *)&(fbin->bin_format)) ) return 1 ;	// read as a 4 byte string
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;		// then endian correct to 4 bytes int
	sprintf(format,"type:%%%lus",sizeof(fbin->bin_type)) ;
	if( read_parameter(fd,format,(void *)&(fbin->bin_type)) ) return 1 ;	// read as a 4 byte string
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;		// then endian correct to 4 bytes int
	config->bin_type = fbin->bin_type ;	// remember this for alvl blocks
	//printf("debug: make_node_fbin: fbin->bin_type=%s\n",strkey(fbin->bin_type)) ;
	return 0 ;
}

int gen_block_fbin(struct node *node, FILE *outfile)
{
	struct block_fbin *fbin = (struct block_fbin *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;
	if( fwrite(&(fbin->bin_format),sizeof(fbin->bin_format),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;
	if( fwrite(&(fbin->bin_type),sizeof(fbin->bin_type),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


int fixup_data_body(struct node *node)
{
	return 0 ;
}

int dump_block_body(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_BODY)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_body(struct node *list, struct config *config, FILE *fd)		// creates a new node for body block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_BODY ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_body(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_gtag
{
	uint32_t gtag ;		// unknown
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_gtag(struct node *node)
{
	if( node->size < sizeof(struct block_gtag) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_gtag)) ;
		return 1 ;
	}
	struct block_gtag *gtag = (struct block_gtag *)node->data ;
	endian_fixup(&(gtag->gtag),sizeof(gtag->gtag)) ;			// fixup the endian order
	return 0 ;
}

int dump_block_gtag(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_gtag) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_gtag)) ;
		return 1 ;
	}
	struct block_gtag *gtag = (struct block_gtag *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_gtag)) ;
	fprintf(outfile,"gtag:%u\n",gtag->gtag) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_gtag(struct node *list, struct config *config, FILE *fd)		// creates a new node for gtag block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_gtag ;
	struct block_gtag *gtag = malloc(sizeof(struct block_gtag)) ;
	if( gtag == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(gtag,0,sizeof(struct block_gtag)) ;
	newnode->data = (unsigned char *)gtag ;
	newnode->size = sizeof(struct block_gtag) ;
	if( read_parameter(fd,"gtag:%u",(void *)&(gtag->gtag)) ) return 1 ;
	return 0 ;
}

int gen_block_gtag(struct node *node, FILE *outfile)
{
	struct block_gtag *gtag = (struct block_gtag *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(gtag->gtag),sizeof(gtag->gtag)) ;
	if( fwrite(&(gtag->gtag),sizeof(gtag->gtag),1,outfile) != 1 ) return 1 ;
	return 0 ;
}

struct block_atag
{
	uint32_t atag ;		// unknown
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_atag(struct node *node)
{
	if( node->size < sizeof(struct block_atag) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_atag)) ;
		return 1 ;
	}
	struct block_atag *atag = (struct block_atag *)node->data ;
	endian_fixup(&(atag->atag),sizeof(atag->atag)) ;			// fixup the endian order
	return 0 ;
}

int dump_block_atag(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_atag) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_atag)) ;
		return 1 ;
	}
	struct block_atag *atag = (struct block_atag *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_atag)) ;
	fprintf(outfile,"atag:%u\n",atag->atag) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_atag(struct node *list, struct config *config, FILE *fd)		// creates a new node for atag block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_atag ;
	struct block_atag *atag = malloc(sizeof(struct block_atag)) ;
	if( atag == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(atag,0,sizeof(struct block_atag)) ;
	newnode->data = (unsigned char *)atag ;
	newnode->size = sizeof(struct block_atag) ;
	if( read_parameter(fd,"atag:%u",(void *)&(atag->atag)) ) return 1 ;
	return 0 ;
}

int gen_block_atag(struct node *node, FILE *outfile)
{
	struct block_atag *atag = (struct block_atag *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(atag->atag),sizeof(atag->atag)) ;
	if( fwrite(&(atag->atag),sizeof(atag->atag),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_indx
{
	uint32_t index ;		// index
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_indx(struct node *node)
{
	if( node->size < sizeof(struct block_indx) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_indx)) ;
		return 1 ;
	}
	struct block_indx *indx = (struct block_indx *)node->data ;
	endian_fixup(&(indx->index),sizeof(indx->index)) ;			// fixup the endian order
	return 0 ;
}

int dump_block_indx(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_indx) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_indx)) ;
		return 1 ;
	}
	struct block_indx *indx = (struct block_indx *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_indx)) ;
	fprintf(outfile,"index:%u\n",indx->index) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_indx(struct node *list, struct config *config, FILE *fd)		// creates a new node for indx block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_indx ;
	struct block_indx *indx = malloc(sizeof(struct block_indx)) ;
	if( indx == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(indx,0,sizeof(struct block_indx)) ;
	newnode->data = (unsigned char *)indx ;
	newnode->size = sizeof(struct block_indx) ;
	if( read_parameter(fd,"index:%u",(void *)&(indx->index)) ) return 1 ;
	return 0 ;
}

int gen_block_indx(struct node *node, FILE *outfile)
{
	struct block_indx *indx = (struct block_indx *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(indx->index),sizeof(indx->index)) ;
	if( fwrite(&(indx->index),sizeof(indx->index),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_scal
{
	double scalar_one ;		// scaling value for I samples
	double scalar_two ;		// scaling value for Q samples
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_scal(struct node *node)
{
	if( node->size < sizeof(struct block_scal) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_scal)) ;
		return 1 ;
	}
	struct block_scal *scal = (struct block_scal *)node->data ;
	endian_fixup(&(scal->scalar_one),sizeof(scal->scalar_one)) ;
	endian_fixup(&(scal->scalar_two),sizeof(scal->scalar_two)) ;
	return 0 ;
}

int dump_block_scal(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_scal) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_scal)) ;
		return 1 ;
	}
	struct block_scal *scal = (struct block_scal *)(node->data) ;
	config->scalar_one = scal->scalar_one ;	// remember this for alvl blocks
	config->scalar_two = scal->scalar_two ;	// remember this for alvl blocks
	fprintf(outfile,"%s\n",strkey(KEY_scal)) ;
	fprintf(outfile,"scalar_one:%.20lf\n",scal->scalar_one) ;
	fprintf(outfile,"scalar_two:%.20lf\n",scal->scalar_two) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_scal(struct node *list, struct config *config, FILE *fd)		// creates a new node for scal block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_scal ;
	struct block_scal *scal = malloc(sizeof(struct block_scal)) ;
	if( scal == NULL )
	{
		printf("Malloc error on data block\n") ;
		return 1 ;
	}
	memset(scal,0,sizeof(struct block_scal)) ;
	newnode->data = (unsigned char *)scal ;
	newnode->size = sizeof(struct block_scal) ;
	if( read_parameter(fd,"scalar_one:%lf",(void *)&(scal->scalar_one)) ) return 1 ;
	if( read_parameter(fd,"scalar_two:%lf",(void *)&(scal->scalar_two)) ) return 1 ;
	config->scalar_one = scal->scalar_one ;	// remember this for alvl blocks
	config->scalar_two = scal->scalar_two ;	// remember this for alvl blocks
	return 0 ;
}

int gen_block_scal(struct node *node, FILE *outfile)
{
	struct block_scal *scal = (struct block_scal *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(scal->scalar_one),sizeof(scal->scalar_one)) ;
	if( fwrite(&(scal->scalar_one),sizeof(scal->scalar_one),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(scal->scalar_two),sizeof(scal->scalar_two)) ;
	if( fwrite(&(scal->scalar_two),sizeof(scal->scalar_two),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_alvl
{
	int16_t isample ;		// I sample
	int16_t qsample ;		// Q sample
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_alvl(struct node *node)
{
	if( node->size < sizeof(struct block_alvl) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_alvl)) ;
		return 1 ;
	}
	struct block_alvl *alvl = (struct block_alvl *)node->data ;
	int nsamples = (node->size)/sizeof(struct block_alvl) ;
	for( int loop = 0 ; loop < nsamples ; loop++, alvl++ )
	{
		endian_fixup(&(alvl->isample),sizeof(alvl->isample)) ;
		endian_fixup(&(alvl->qsample),sizeof(alvl->qsample)) ;
	}
	return 0 ;
}

int dump_block_alvl(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_alvl) )
	{
		printf("Block '%s' is truncated\n",strkey(KEY_alvl)) ;
		return 1 ;
	}
	double factor = 1.0L ;
	switch ( (uint32_t )config->bin_type )
	{
		case (uint32_t )BINTYPE_FLT4:
			factor = (double )1 ;
		break ;
		case (uint32_t )BINTYPE_FIX2:
			factor = (double )0x7FFF ;
		break ;
		case (uint32_t )BINTYPE_FIX3:
			factor = (double )0x7FFFFF ;
		break ;
		case (uint32_t )BINTYPE_FIX4:
			factor = (double )0x7FFFFFFF ;
		break ;
		default:
			printf("Unknown bin_type '%s' (%x) at index %d\n",strkey(config->bin_type),config->bin_type,config->index) ;
			return 1;
	}
	fprintf(outfile,"%s\n",strkey(KEY_alvl)) ;
	struct block_alvl *alvl = (struct block_alvl *)(node->data) ;
	int nsamples = (node->size)/sizeof(struct block_alvl) ;
	for( int loop = 0 ; loop < nsamples ; loop++, alvl++ )
	{
		int isample = alvl->isample ;
		int qsample = alvl->qsample ;
		double scaled_i = (double )isample/factor*config->scalar_one ;	// precedence order?
		double scaled_q = (double )qsample/factor*config->scalar_two ;
		fprintf(outfile,"i:%.20lf\n",scaled_i) ;
		fprintf(outfile,"q:%.20lf\n",scaled_q) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int read_alvl_samples(struct block_alvl *, int, struct config *, FILE *) ;

int make_node_alvl(struct node *list, struct config *config, FILE *fd)		// creates a new node for alvl block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_alvl ;
	int alvl_lines = count_alvl_lines(fd) ;		// count lines, 2 lines per sample (i and q), use this to malloc space for the entire block
	if( alvl_lines <= 0 )
	{
		printf("Error counting lines in '%s' block\n",strkey(KEY_alvl)) ;
		return 1 ;
	}
	if( alvl_lines % 2 == 1 )
	{
		printf("Cannot deal with an odd number of lines (%d) reading '%s' block\n",alvl_lines,strkey(KEY_alvl)) ;
		return 1 ;
	}
	//printf("debug: make_node_alvl: alvl_lines=%u\n",alvl_lines) ;
	int alvl_samples = alvl_lines/2 ;	// a sample is a pair of I,Q values
	size_t alvl_size = alvl_samples * sizeof(struct block_alvl) ;
	//printf("debug: make_node_alvl: alvl_samples=%d alvl_size=%zu\n",alvl_samples,alvl_size) ;
	struct block_alvl *alvl_data = malloc(alvl_size) ;
	if( alvl_data == NULL )
	{
		printf("Malloc error on '%s' data block\n",strkey(KEY_alvl)) ;
		return 1 ;
	}
	memset(alvl_data,0,alvl_size) ;
	newnode->data = (unsigned char *)alvl_data ;
	newnode->size = alvl_size ;
	if( read_alvl_samples(alvl_data,alvl_samples,config,fd) )	// read lines of i, q values, convert them into ints and store them in alvl_data
	{
		printf("Error reading '%s' block\n",strkey(KEY_alvl)) ;
		return 1 ;
	}
	return 0 ;
}

int count_alvl_lines(FILE *fd)		// count how many lines of i,q data there are before we hit a blank line (end of block)
{
	char line[SIZE_LINE] ;
	int count = 0 ;
	unsigned long alvl_start = ftell(fd) ;
	while( fgets(line,SIZE_LINE,fd) != NULL )
	{
		chomp(line,SIZE_LINE) ;				// remove newline
		if( strlen(line) < 1 )
			break ;
		count++ ;
	}
	fseek(fd,alvl_start,SEEK_SET) ;
	return count ;
}

int read_alvl_samples(struct block_alvl *alvl_data, int alvl_samples, struct config *config, FILE *fd)
{
	char line1[SIZE_LINE] ;
	char line2[SIZE_LINE] ;
	double factor = 1.0L ;
	unsigned long alvl_start = ftell(fd) ;
	switch ( (uint32_t )config->bin_type )
	{
		case (uint32_t )BINTYPE_FLT4:
			factor = (double )1 ;
		break ;
		case (uint32_t )BINTYPE_FIX2:
			factor = (double )0x7FFF ;
		break ;
		case (uint32_t )BINTYPE_FIX3:
			factor = (double )0x7FFFFF ;
		break ;
		case (uint32_t )BINTYPE_FIX4:
			factor = (double )0x7FFFFFFF ;
		break ;
		default:
			printf("Unknown bin_type '%s' (%x) at index %d\n",strkey(config->bin_type),config->bin_type,config->index) ;
			return 1;
	}
	for( int sample_count = 0 ; sample_count < alvl_samples ; sample_count++, alvl_data++ )
	{
		if( fgets(line1,SIZE_LINE,fd) == NULL ) return 1 ;
		if( fgets(line2,SIZE_LINE,fd) == NULL ) return 1 ;
		chomp(line1,SIZE_LINE) ; chomp(line2,SIZE_LINE) ;
		if( strlen(line1) == 0 ) return 1 ;
		if( strlen(line2) == 0 ) return 1 ;
		double i ; double q ;
		int convert_count = sscanf(line1,"i:%lf",&i) ;
		if( convert_count != 1 )
		{
			printf("Failed to read 'i' value %d from line %s\n",sample_count,line1) ;
			return 1 ;
		}
		convert_count = sscanf(line2,"q:%lf",&q) ;
		if( convert_count != 1 )
		{
			printf("Failed to read 'q' value %d from line %s\n",sample_count,line2) ;
			return 1 ;
		}
		double i_scaled = (i/config->scalar_one)*factor ;
		double q_scaled = (q/config->scalar_two)*factor ;
		int i_int = round(i_scaled) ;
		int q_int = round(q_scaled) ;
		alvl_data->isample = i_int ;
		alvl_data->qsample = q_int ;
		//if( sample_count == 0 ) printf("debug: read_alvl_samples: double i=%lf q=%lf, scalar_one=%lf scalar_two=%lf, factor=%lf int i=%d q=%d\n",i,q,config->scalar_one,config->scalar_two,factor,alvl_data->isample,alvl_data->qsample) ;
	}
	fseek(fd,alvl_start,SEEK_SET) ;
	return 0 ;
}

int gen_block_alvl(struct node *node, FILE *outfile)
{
	struct block_alvl *alvl = (struct block_alvl *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	int actual_size = node->size ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	int sample_count = actual_size/sizeof(struct block_alvl) ;
	//printf("debug: gen_block_alvl: actual size %d, sample_count %d\n",actual_size,sample_count) ;
	for( int count = 0 ; count < sample_count ; count++, alvl++ )
	{
		//if( count == 0 ) printf("debug: gen_block_alvl: sample 0 i=%d q=%d\n",alvl->isample,alvl->qsample) ;
		endian_fixup(&(alvl->isample),sizeof(alvl->isample)) ;
		if( fwrite(&(alvl->isample),sizeof(alvl->isample),1,outfile) != 1 ) return 1 ;
		endian_fixup(&(alvl->qsample),sizeof(alvl->qsample)) ;
		if( fwrite(&(alvl->qsample),sizeof(alvl->qsample),1,outfile) != 1 ) return 1 ;
	}
	return 0 ;
}


int fixup_data_end(struct node *node)
{
	return 0 ;
}

int dump_block_end(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_END)) ;
	return 0 ;
}

int make_node_end(struct node *list, struct config *config, FILE *fd)		// creates a new node for end block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		printf("Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_END ;
	newnode->size = 0 ;
	newnode->data = NULL ;
	return 0 ;
}

int gen_block_end(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


void free_all_nodes(struct node *list)
{
	while( list != NULL )
	{
		struct node *next = list->next ;
		free(list) ;
		list = next ;
	}
}

void free_all_nodes_and_data(struct node *list)
{
	while( list != NULL )
	{
		if( list->data != NULL )
			free(list->data) ;
		struct node *next = list->next ;
		free(list) ;
		list = next ;
	}
}

//END
