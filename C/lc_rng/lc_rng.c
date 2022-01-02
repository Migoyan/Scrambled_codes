/**
 * @file lc_rng.c
 * @brief lc_rng generator
 * 
 */
#include <stdio.h>
#include <math.h>

#define N 100000
#define x_0 12

const long k = 16807;
const long l = 0;
const long m = 2147483647;

/**
 * @brief Linear congruent generator
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

int main()
{
	long x[N];
	x[0] = x_0;
	FILE *file = fopen("./result/rng.dat", "w");
	for (long i = 0; i < N - 1; i++)
	{
		x[i + 1] = lc_rng(x[i], k, l, m);
		fprintf(file, "%ld;%ld\n", x[i], x[i + 1]);
	}
	fclose(file);

	return 0;
}
