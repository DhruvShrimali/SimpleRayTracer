#include <cstdlib>
#include <cmath>
#include "random.h"


float randomFloat(){
    // Convert to [-0.5, 0.5] range
    return (float)(((rand() % 1000)-500)/1000.0f);
}

float randomProbability(){
    return (float)((rand() % 10000)/10000.0f);
}

float guassianRandomFloat() {
    // Convert to [0, 1] range
    double u1 = (rand() % 1000)/1000.0f;
    double u2 = (rand() % 1000)/1000.0f;

    // Box-Muller transform
    double z0 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);  // Gaussian distributed number

    // Convert to [-0.5, 0.5] range (mostly)
    return z0/8.0f;
}