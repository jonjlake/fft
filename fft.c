#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "fft.h"

double complex_magnitude_squared(ComplexNo num_in)
{
	return (pow(num_in.re, 2.0) + pow(num_in.im, 2.0));
}

double complex_magnitude(ComplexNo num_in)
{
	return sqrt(complex_magnitude_squared(num_in));
}

ComplexNo *fourier_transform(double *array_in, int num_samples)
{
	ComplexNo *array_out = (ComplexNo*)calloc(num_samples, sizeof(*array_out));
	int i, j;

	for (i = 0; i < num_samples; i++)
	{

		for (j = 0; j < num_samples; j++)
		{

			array_out[i].re += array_in[j] * cos(2.0 * M_PI * j * i / (double)num_samples);
			array_out[i].im -= array_in[j] * sin(2.0 * M_PI * j * i / (double)num_samples);
		}
	}

	return array_out;
}

double *fourier_power(ComplexNo *fourier_array, int num_samples)
{
	double *output_array = (double *)calloc(num_samples, sizeof(*output_array));
	int i;

	for (i = 0; i < num_samples; i++)
	{
		output_array[i] = complex_magnitude_squared(fourier_array[i]);
	}

	return output_array;
}

void fft_1(void)
{
	double sample_array[5] = { 0, 1, 0, -1, 0};
	double sample_freq = 0.2; // Hz
	int num_samples = 5;

	double fourier_frequencies[5] = { 0 };
	ComplexNo *fourier_transform_array = NULL;

	int i;

	for (i = 0; i < num_samples; i++)
	{
		fourier_frequencies[i] = (double) i / (M_PI * (double)num_samples);
	}

	fourier_transform_array = fourier_transform(sample_array, num_samples);

	for (i = 0; i < num_samples; i++)
	{
		printf("Frequency %f: Re = %f, Im = %f\n", fourier_frequencies[i], fourier_transform_array[i].re, fourier_transform_array[i].im);
	}

	free(fourier_transform_array);
}

double *generate_sine(double A, double f, double dt, double ph, int num_samples)
{
	double *output_array = (double *)calloc(num_samples, sizeof(*output_array));
	int i;

	for (i = 0; i < num_samples; i++)
	{
		output_array[i] = A * sin(2.0 * M_PI * f * i * dt + ph);
	}

	return output_array;
}

double *generate_fourier_frequencies(double dt, int num_samples)
{
	double *output_array = (double *)calloc(num_samples, sizeof(*output_array));
	int i;

	for (i = 0; i < num_samples; i++)
	{
		output_array[i] = (double)i / (dt * num_samples);
	}

	return output_array;
}

void print_ft(double *ft_frequencies, ComplexNo *ft, int num_samples)
{
	int i;

	for (i = 0; i < num_samples; i++)
	{
		printf("Frequency %f: Re = %f, Im = %f\n", ft_frequencies[i], ft[i].re, ft[i].im);
	}
}

void print_ft_to_csv(char *filename, double dt, double *sample_array, double *sample_frequencies, ComplexNo *sample_ft, double *sample_power, int num_samples)
{
	FILE *fp = fopen(filename, "w+");
	int i;

	fprintf(fp, "Time, Value, Frequency, FT Re, FT Im, FT Mag^2\n");

	for (i = 0; i < num_samples; i++)
	{
		fprintf(fp, "%f, %f, %f, %f, %f, %f\n", dt * (double)i, sample_array[i], sample_frequencies[i], 
				sample_ft[i].re, sample_ft[i].im, sample_power[i]);
	}	
}

void fft_2(void)
{
	double fsamp = 1.0;
	double f = 0.2;
	double dt = 1;
	int num_samples = 50;
	double *sample_array = generate_sine(1, f, dt, 0, num_samples);
	double *sample_frequencies = generate_fourier_frequencies(dt, num_samples);
	ComplexNo *sample_ft = fourier_transform(sample_array, num_samples);
	double *sample_power = fourier_power(sample_ft, num_samples);

	print_ft(sample_frequencies, sample_ft, num_samples);
	print_ft_to_csv("fft_2.csv", dt, sample_array, sample_frequencies, sample_ft, sample_power, num_samples);
	free(sample_array);
	free(sample_frequencies);
	free(sample_ft);
	free(sample_power);
}

void do_sine_ft(char *filename, double A, double fsamp, double f, double dt, double ph, int num_samples) 
{
	double *sample_array = generate_sine(A, f, dt, ph, num_samples);
	double *sample_frequencies = generate_fourier_frequencies(dt, num_samples);
	ComplexNo *sample_ft = fourier_transform(sample_array, num_samples);
	double *sample_power = fourier_power(sample_ft, num_samples);

	print_ft(sample_frequencies, sample_ft, num_samples);
	print_ft_to_csv(filename, dt, sample_array, sample_frequencies, sample_ft, sample_power, num_samples);
	free(sample_array);
	free(sample_frequencies);
	free(sample_ft);
	free(sample_power);
}

void fft_3(void)
{
	do_sine_ft("fft_3.csv", 1, 1, 0.2, 1, M_PI / 2.0, 50);
}

void fft_4(void)
{
	do_sine_ft("fft_4.csv", 1, 1, 0.2123, 1, 0, 50);
}

int main()
{
	fft_1();
	fft_2();
	fft_3();
	fft_4();

	return 0;
}

