/**
 *  @file testfloat_CompDecomp.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sz.h"
#include "rw.h"
#include "zc.h"

int main(int argc, char * argv[])
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char outDir[640], oriFilePath[640], outputFilePath[640];
    char *cfgFile, *zcFile, *solName, *varName, *errBoundMode;
    double absErrBound;
    int errboundmode;
    if(argc < 2)
    {
	printf("Test case: testfloat_CompDecomp2 [config_file] [srcFilePath] [dimension sizes...]\n");
	printf("Example: testfloat_CompDecomp2 sz.config testdata/x86/testfloat_8_8_128.dat 8 8 128\n");
	exit(0);
    }
   
    cfgFile=argv[1];

    sprintf(oriFilePath, "%s", argv[2]);
    if(argc>=4)
	r1 = atoi(argv[3]); //8
    if(argc>=5)
	r2 = atoi(argv[4]); //8
    if(argc>=6)
	r3 = atoi(argv[5]); //128
    if(argc>=7)
        r4 = atoi(argv[6]);
    if(argc>=8)
        r5 = atoi(argv[7]);
   
    printf("cfgFile=%s\n", cfgFile); 
    SZ_Init(cfgFile);
   
    //printf("zcFile=%s\n", zcFile);
    //ZC_Init(zcFile);
 
    sprintf(outputFilePath, "%s.sz", oriFilePath);
  
    size_t nbEle; 
    int status = SZ_SCES;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
   
    size_t outSize; 
    printf("1\n");
    unsigned char *bytes = SZ_compress(SZ_FLOAT, data, &outSize, r5, r4, r3, r2, r1);
    float *decData = SZ_decompress(SZ_FLOAT, bytes, outSize, r5, r4, r3, r2, r1);
    free(bytes);
    free(decData);
    free(data);

    data = readFloatData(oriFilePath, &nbEle, &status);
    printf("2\n");
    bytes = SZ_compress(SZ_FLOAT, data, &outSize, r5, r4, r3, r2, r1);
    decData = SZ_decompress(SZ_FLOAT, bytes, outSize, r5, r4, r3, r2, r1);
    free(bytes);
    free(decData);
    free(data);
        
    printf("3\n");
    
    printf("done\n");
    
    SZ_Finalize();
    //ZC_Finalize();
    return 0;
}
