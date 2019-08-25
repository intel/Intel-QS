#include <stdio.h> 
#include <stdlib.h>
#include "hdf5.h"

#define BINFILE "testfloat_8_8_128.dat"
#define HDF5FILE "testfloat_8_8_128.h5"
int main() {

	hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
	hsize_t     dims[3];
	herr_t      status;

	int nbEle = 8*8*128;
	FILE *f;
	f = fopen(BINFILE, "rb");
	float *data = (float*)malloc(nbEle*sizeof(float));
	fread(data, sizeof(float), nbEle, f);
	fclose(f);

	/* Create a new file using default properties. */
	file_id = H5Fcreate(HDF5FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Create the data space for the dataset. */
	dims[0] = 128; 
	dims[1] = 8; 
	dims[2] = 8;
	dataspace_id = H5Screate_simple(3, dims, NULL);

	/* Create the dataset. */
	dataset_id = H5Dcreate2(file_id, "/testfloat", H5T_IEEE_F32LE, dataspace_id, 
						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
					 data);

	/* End access to the dataset and release resources used by it. */
	status = H5Dclose(dataset_id);

	/* Terminate access to the data space. */ 
	status = H5Sclose(dataspace_id);

	/* Close the file. */
	status = H5Fclose(file_id);
}
