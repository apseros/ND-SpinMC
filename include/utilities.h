#ifndef UTILITIES_H
#define UTILITIES_H

#include <random>

extern std::mt19937 rng;

template<typename T>
int sgn(T val) { return (T(0.0) < val) - (val < T(0.0)); }

void random_sampling(const double *, double *, int);

double mean(const double *, int);

double stdev(const double *, int);

double bootstrap(const double *, int, double (*)(const double *, int), int);

int int_pow(int, int);

double correlation(const double *, int, int);

void autocorrelation(const double *, double *, int);

#endif
