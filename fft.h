/* This file could use some structs to package data */

#define FREE_AND_NULL(x) if (NULL != x) {free(x); x = NULL;}

typedef struct ComplexNo{
	double re;
	double im;
} ComplexNo;

typedef struct FourierData{
	double *sample_array;
	double *sample_frequencies;
	ComplexNo *sample_ft;
	ComplexNo *sample_ift;
	double *sample_power;
	double dt;
	int num_samples;
} FourierData;

ComplexNo *fourier_transform(double *array_in, int num_samples);

void calculate_ft_arrays(FourierData *p_fourier_data);
void calculate_ift_arrays(FourierData *p_fourier_data);

void destroy_fourier_data_arrays(FourierData *p_fourier_data);

