/**
 *  @file H5Z_SZ.c
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief SZ filter for HDF5
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "H5Z_SZ.h"
#include "H5PLextern.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
        totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


//sz_params* conf_params = NULL;

int load_conffile_flag = 0; //0 means 'not yet', 1 means 'already loaded'
char cfgFile[256] = "sz.config"; 

const H5Z_class2_t H5Z_SZ[1] = {{
	H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
	(H5Z_filter_t)H5Z_FILTER_SZ, /* Filter id number */
	1,              /* encoder_present flag (set to true) */
	1,              /* decoder_present flag (set to true) */
	"SZ compressor/decompressor for floating-point data.", /* Filter name for debugging */
	NULL,                          /* The "can apply" callback */
	H5Z_sz_set_local,                          /* The "set local" callback */
	(H5Z_func_t)H5Z_filter_sz,   /* The actual filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_SZ;}

int H5Z_SZ_Init(char* cfgFile) 
{ 
	herr_t ret;
	//printf("start in H5Z_SZ_Init, load_conffile_flag = %d\n", load_conffile_flag);
	if(load_conffile_flag==0)
	{
		load_conffile_flag = 1;
		int status = SZ_Init(cfgFile);
		//printf("cfgFile=%s\n", cfgFile);
		//printf("szMode=%d, errorBoundMode=%d, relBoundRatio=%f\n", szMode, errorBoundMode, relBoundRatio);
		if(status == SZ_NSCS)
			return SZ_NSCS;
		else
			return SZ_SCES;		
	}

	ret = H5Zregister(H5Z_SZ); 
	if(ret < 0)
		return SZ_NSCS;
	else
		return SZ_SCES;
}

int H5Z_SZ_Init_Params(sz_params *params) 
{ 
	herr_t ret = H5Zregister(H5Z_SZ); 
	int status = SZ_Init_Params(params);
	if(status == SZ_NSCS || ret < 0)
		return SZ_NSCS;
	else
		return SZ_SCES;
}

sz_params* H5Z_SZ_Init_Default()
{
	herr_t ret = H5Zregister(H5Z_SZ);	
	
	sz_params* conf_params = (sz_params *)malloc(sizeof(sz_params));
	confparams_cpr->quantization_intervals = 0;
	confparams_cpr->max_quant_intervals = 65536;
    dataEndianType = LITTLE_ENDIAN_DATA;
    confparams_cpr->sol_ID = SZ;
    //confparams_cpr->layers = 1;
    confparams_cpr->sampleDistance = 100;
    confparams_cpr->predThreshold = 0.99;
    //confparams_cpr->offset = 0;
    confparams_cpr->szMode = SZ_BEST_COMPRESSION;
    confparams_cpr->gzipMode = 1; //best speed
    confparams_cpr->errorBoundMode = REL; //details about errorBoundMode can be found in sz.config
    confparams_cpr->absErrBound = 1E-4;
    confparams_cpr->relBoundRatio = 1E-3;
    confparams_cpr->pw_relBoundRatio = 1E-4;
    confparams_cpr->segment_size = 32;
    confparams_cpr->pwr_type = SZ_PWR_AVG_TYPE;	
	
	int status = SZ_Init_Params(conf_params);
	if(status == SZ_NSCS || ret < 0)
		return NULL;
	else
		return confparams_cpr;
}

int H5Z_SZ_Finalize()
{
	//SZ_Finalize();
	herr_t ret = H5Zunregister(H5Z_FILTER_SZ);
	if (ret < 0) return -1;
	return 0;
}

/**
 * to be used in decompression and compression, inside the H5Z_filter_sz().
 * */
void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1)
{
	//printf("cd_nelmts=%zu\n", cd_nelmts); 
	assert(cd_nelmts >= 4);
	unsigned char bytes[8];	
	*dimSize = cd_values[0];
	*dataType = cd_values[1];

	switch(*dimSize)
	{
	case 1:
		intToBytes_bigEndian(bytes, cd_values[2]);
		intToBytes_bigEndian(&bytes[4], cd_values[3]);
		if(sizeof(size_t)==4)
			*r1 = (unsigned int)bytesToLong_bigEndian(bytes);
		else
			*r1 = (unsigned long)bytesToLong_bigEndian(bytes);
		*r2 = *r3 = *r4 = *r5 = 0;
		break;
	case 2:
		*r3 = *r4 = *r5 = 0;
		*r2 = cd_values[3];
		*r1 = cd_values[2];
		break;
	case 3:
		*r4 = *r5 = 0;
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];
		break;
	case 4:
		*r5 = 0;
		*r4 = cd_values[5];
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];	
		break;
	default: 
		*r5 = cd_values[6];
		*r4 = cd_values[5];
		*r3 = cd_values[4];
		*r2 = cd_values[3];
		*r1 = cd_values[2];		
	}
	//printf("SZ_cdArrayToMetaData: r5=%zu, r4=%zu, r3=%zu, r2=%zu, r1=%zu\n", *r5, *r4, *r3, *r2, *r1);	
}

/**
 * to be used in compression, and to be called outside H5Z_filter_sz().
 * */
void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char bytes[8] = {0};
	unsigned long size;
	*cd_values = (unsigned int*)malloc(sizeof(unsigned int)*7);
	int dim = computeDimension(r5, r4, r3, r2, r1);
	(*cd_values)[0] = dim;
	(*cd_values)[1] = dataType;	//0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....
	switch(dim)
	{
	case 1:
		size = (unsigned long)r1;
		longToBytes_bigEndian(bytes, size);
		(*cd_values)[2] = bytesToInt_bigEndian(bytes);
		(*cd_values)[3] = bytesToInt_bigEndian(&bytes[4]);	
		*cd_nelmts = 4;
		break;
	case 2:
		(*cd_values)[2] = (unsigned int) r2;
		(*cd_values)[3] = (unsigned int) r1;	
		*cd_nelmts = 4;
		break;
	case 3:
		(*cd_values)[2] = (unsigned int) r3;
		(*cd_values)[3] = (unsigned int) r2;
		(*cd_values)[4] = (unsigned int) r1;	
		*cd_nelmts = 5;
		break;
	case 4:
		(*cd_values)[2] = (unsigned int) r4;	
		(*cd_values)[3] = (unsigned int) r3;
		(*cd_values)[4] = (unsigned int) r2;
		(*cd_values)[5] = (unsigned int) r1;	
		*cd_nelmts = 6;
		break;
	default:
		(*cd_values)[2] = (unsigned int) r5;		
		(*cd_values)[3] = (unsigned int) r4;	
		(*cd_values)[4] = (unsigned int) r3;
		(*cd_values)[5] = (unsigned int) r2;
		(*cd_values)[6] = (unsigned int) r1;
		*cd_nelmts = 7;	
	}
}

static herr_t H5Z_sz_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id)
{
	//printf("start in H5Z_sz_set_local\n");
	size_t r5=0,r4=0,r3=0,r2=0,r1=0, dsize;
	static char const *_funcname_ = "H5Z_sz_set_local";
	int i, ndims, ndims_used = 0;	
	hsize_t dims[H5S_MAX_RANK], dims_used[5] = {0,0,0,0,0};	
	herr_t retval = 0;
	H5T_class_t dclass;
	H5T_sign_t dsign;
	unsigned int flags = 0;
	//conf_params = H5Z_SZ_Init_Default();
	H5Z_SZ_Init(cfgFile);
	
	int dataType = SZ_FLOAT;
	
	if (0 > (dclass = H5Tget_class(type_id)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

	if (0 == (dsize = H5Tget_size(type_id)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

	if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims, 0)))
		H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
		
	for (i = 0; i < ndims; i++)
	{
		if (dims[i] <= 1) continue;
		dims_used[ndims_used] = dims[i];
		ndims_used++;
	}
	
	//printf("dclass=%d, H5T_FLOAT=%d, H5T_INTEGER=%d\n", dclass, H5T_FLOAT, H5T_INTEGER);
	if (dclass == H5T_FLOAT)
		dataType = dsize==4? SZ_FLOAT: SZ_DOUBLE;
	else if(dclass == H5T_INTEGER)
	{
		if (0 > (dsign = H5Tget_sign(type_id)))
			H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");		
		if(dsign == H5T_SGN_NONE) //unsigned
		{
			switch(dsize)
			{
			case 1:
				dataType = SZ_UINT8;
				break;
			case 2:
				dataType = SZ_UINT16;
				break;
			case 4:
				dataType = SZ_UINT32;
				break;
			case 8:
				dataType = SZ_UINT64;
				break;
			}
		}
		else
		{
			switch(dsize)
			{
			case 1:
				dataType = SZ_INT8;
				break;
			case 2:
				dataType = SZ_INT16;
				break;
			case 4:
				dataType = SZ_INT32;
				break;
			case 8:
				dataType = SZ_INT64;
				break;
			}			
		}
	}
	else
	{
		//printf("Error: dclass...\n");
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
	}
	
	
	switch(ndims_used)
	{
	case 1: 
		r1 = dims_used[0];
		break;
	case 2:
		r1 = dims_used[0];
		r2 = dims_used[1];
		break;
	case 3:
		r1 = dims_used[0];
		r2 = dims_used[1];
		r3 = dims_used[2];		
		break;
	case 4:
		r1 = dims_used[0];
		r2 = dims_used[1];
		r3 = dims_used[2];	
		r4 = dims_used[3];
		break;
	default: 
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "requires chunks w/1,2,3 or 4 non-unity dims");
	}
	
	size_t cd_nelmts = 0;
	unsigned int mem_cd_values[7]; 
	unsigned int* cd_values;

	if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_SZ, &flags, &cd_nelmts, mem_cd_values, 0, NULL, NULL))
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current SZ cd_values");

	SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1);
	
	/* Now, update cd_values for the filter */
	if (0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_SZ, flags, cd_nelmts, cd_values))
		H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");	
		
	retval = 1;
done:
	return retval;
}


static size_t H5Z_filter_sz(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf)
{
//	printf("start in H5Z_filter_sz\n");
	//H5Z_SZ_Init_Default();
	
	size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
	int dimSize = 0, dataType = 0;

	if(cd_nelmts==0) //this is special data such as string, which should not be treated as values.
		return nbytes;
	
	SZ_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1);
	
/*	int i=0;
	for(i=0;i<cd_nelmts;i++)
		printf("cd_values[%d]=%u\n", i, cd_values[i]);
	printf("dimSize=%d, r1=%u, r2=%u, r3=%u, r4=%u, r5=%u\n", dimSize, r1, r2, r3, r4, r5);*/
	size_t nbEle = computeDataLength(r5, r4, r3, r2, r1); 
	
	if (flags & H5Z_FLAG_REVERSE) 
	{ 
		//cost_start();
		/* decompress data */
		if(dataType == SZ_FLOAT)//==0
		{
			float* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);

			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(float);
		}
		else if(dataType == SZ_DOUBLE)//==1
		{
			double* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(double);			
		}
		else if(dataType == SZ_INT8)
		{
			char* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(char);
		}
		else if(dataType == SZ_UINT8)
		{
			unsigned char* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned char);
		}
		else if(dataType == SZ_INT16)
		{
			short* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(short);
		}
		else if(dataType == SZ_UINT16)
		{
			unsigned short* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned short);
		}
		else if(dataType == SZ_INT32)
		{
			int* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(int);
		}
		else if(dataType == SZ_UINT32)
		{
			unsigned int* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned int);
		}
		else if(dataType == SZ_INT64)
		{
			long* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(long);
		}
		else if(dataType == SZ_UINT64)
		{
			unsigned long* data = SZ_decompress(dataType, *buf, nbytes, r5, r4, r3, r2, r1);
										
			free(*buf);
			*buf = data;
			*buf_size = nbEle*sizeof(unsigned long);
		}
		else
		{
			printf("Decompression error: unknown data type: %d\n", dataType);
			exit(0);
		}
		//cost_end();
		//printf("decompression time = %lf, decompression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
	}
	else
	{
		size_t outSize = 0;
	
		//printf("r5=%d, r4=%d, r3=%d, r2=%d, r1=%d, dataType=%d\n", r5, r4, r3, r2, r1, dataType);
		//cost_start();
		if(dataType == SZ_FLOAT)//==0
		{
			float* data = (float*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_DOUBLE)//==1
		{
			double* data = (double*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_INT8)
		{
			char* data = (char*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_UINT8)
		{
			unsigned char* data = (unsigned char*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_INT16)
		{
			short* data = (short*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_UINT16)
		{
			unsigned short* data = (unsigned short*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_INT32)
		{
			int* data = (int*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_UINT32)
		{
			unsigned int* data = (unsigned int*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_INT64)
		{
			long* data = (long*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else if(dataType == SZ_UINT64)
		{
			unsigned long* data = (unsigned long*)(*buf);
			unsigned char *bytes = SZ_compress(dataType, data, &outSize, r5, r4, r3, r2, r1);
			free(*buf);
			*buf = bytes;
			*buf_size = outSize;
		}
		else 
		{
			printf("Compression error: unknown data type: %d\n", dataType);
			exit(0);
		}
		//cost_end();
		//printf("compression time = %lf, compression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
	}
	
	H5Z_SZ_Finalize();
	return *buf_size;
}

void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	switch(dim)
	{
	case 1: 
		dims[0] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
			chunk[0] = r1;
		else
			chunk[0] = 2147483648;//2^31
		break;
	case 2:
		dims[0] = r2;
		dims[1] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r2;
			chunk[1] = r1;
		}
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}	
		break;
	case 3:
		dims[0] = r3;
		dims[1] = r2;
		dims[2] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r3;
			chunk[1] = r2;
			chunk[2] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
		break;
	case 4:
		dims[0] = r4;
		dims[1] = r3;
		dims[2] = r2;
		dims[3] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r4;
			chunk[1] = r3;
			chunk[2] = r2;
			chunk[3] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
		break;
	default:
		dims[0] = r5;
		dims[1] = r4;
		dims[2] = r3;
		dims[3] = r2;
		dims[4] = r1;
		if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
		{
			chunk[0] = r5;
			chunk[1] = r4;
			chunk[2] = r3;
			chunk[3] = r2;
			chunk[4] = r1;
		}		
		else
		{
			printf("Error: size is too big!\n");
			exit(0);
		}
	}
}


