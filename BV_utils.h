#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define PARI_OLD_NAMES
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <pari/pari.h>
//#include <pari.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

struct timeval tv;

double Uniform(void) {
    return ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
}

double Normal(void) {
    return sqrt(-log(Uniform())*2.0) * sin(2.0 * M_PI * Uniform());
}

double Gauss(double mu, double sigma) {
    double z = sqrt(-2.0 * log(Uniform())) * sin(2.0 * M_PI * Uniform());
    return mu + sigma*z;
}

GEN Sample(int n, double sigma) {
    GEN ret = cgetg(n + 1, t_VEC);
    double z;
    int i;
    
    for (i = 1; i <= n; i++) {
        z = Gauss(0, sigma);
        z = fabs(round(z)); /*absolute value of Gaussian distribution */
        ret[i] = (long) stoi((long) z);
    }
    
    return ret;
}

GEN generate_random(int bit_length){
    gettimeofday(&tv, NULL);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    GEN r = randomi(gshift(gen_1, bit_length));
    return r;
}
