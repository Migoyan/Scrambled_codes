/**
 * @file autocorrelation.c
 * @brief program that compute autocorrelation function
 * 
 */
#include <stdio.h>
#include <math.h>

#define N 100
#define x_0 12

const long k = 899;
const long l = 0;
const long m = 32768;

/**
 * @brief linear congruential random number generator
 * 
 * @param x 
 * @param k 
 * @param l 
 * @param m 
 * @return long 
 */
static inline long lc_rng(long x, long k, long l, long m)
{
	return (k * x + l) % m;
}

double correlation_rate(long ksi[])
{
	long long ksi_m, ksi_m_1, ksi_0, ksi_1, ksi_produit, ksi_temp, ksi_temp_1;

	for (long i = 0; i < N - 1; i++)
	{
		ksi_m += ksi[i];
		ksi_m_1 += ksi[i+1];
	}
	ksi_m /= N-1;
	ksi_m_1 /= N-1;

	for (long i = 0; i < N - 1; i++)
	{
		ksi_temp = ksi[i] - ksi_m;
		ksi_temp_1 = ksi[i + 1] - ksi_m_1;
		ksi_produit += ksi_temp_1 * ksi_temp;
		ksi_0 += ksi_temp * ksi_temp;
		ksi_1 += ksi_temp_1 * ksi_temp_1;
	}
	return (double) ksi_produit / sqrt((double)ksi_0 * (double)ksi_1);
}

int main()
{
	long x[N];
	x[0] = x_0;
	//FILE *file = fopen("./result/rng.dat", "w");
	for (long i = 0; i < N - 1; i++)
	{
		x[i + 1] = lc_rng(x[i], k, l, m);
	//	fprintf(file, "%ld;%ld\n", x[i], x[i + 1]);
	}
	//fclose(file);

	printf("%0.15lf\n", correlation_rate(x));

	return 0;
}