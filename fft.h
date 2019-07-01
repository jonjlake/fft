/* This file could use some structs to package data */

typedef struct ComplexNo{
	double re;
	double im;
} ComplexNo;

ComplexNo *fourier_transform(double *array_in, int num_samples);


