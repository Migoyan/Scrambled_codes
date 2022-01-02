/**
 * @file lc_rng.c
 * @brief lc_rng generator
 * 
 */
#include <stdio.h>
#include <math.h>

#define k 16807
#define l 0
#define m 2147483647

int main ()
{
    const long x_0 = 12;
    long x = x_0;
    unsigned long i = 0;
    
    // generator
    do
    {
        x = (k * x + l)%m;
        i++;
    }
    while (x_0 != x);
    
    printf("p = %ld\n", i);
    return 0;
}
