#include <stdio.h>
#include <stdlib.h>

// Generate a uniform random number strictly between 0 and 1 (avoids endpoints in initial positions)
// Note: this is NOT the best random number out there, but it's quick and simple
double generate_uniform_random(void)
{
  double r;
  r = ((double)rand() + 0.5) / ((double)RAND_MAX + 1.0);
  return r;
}

// Generate a pseudo-random number from a near-Gaussian distribution with zero average and unit variance.
// Central limit theorem: a sum of 12 uniform random numbers on [0,1] has average 6 and variance 1,
// and is close to normally distributed. (The tails are cut off at +/- 6 standard
// deviations, which is of no consequence for initializing velocities.)
double gauss(void)
{
  double sum = -6.0;
  for (int i = 0; i < 12; i++)
  {
    sum += generate_uniform_random();
  }
  return sum;
}
