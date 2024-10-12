#include <cmath>
#include <iostream>
#include "utilities.h"

std::random_device rd;  // seed
std::mt19937 rng{rd()};


void random_sampling(const double *data, double *result, int array_size) {
    std::uniform_int_distribution<int> distrib{0, array_size - 1};
    for (int i = 0; i < array_size; i++) {
        result[i] = data[distrib(rng)];
    }
}

double bootstrap(const double *data, int array_size, double (*operation)(const double *, int), int bins) {
    auto *bin_stat = new double[bins]();
    for (int i = 0; i < bins; i++) {
        // draw random sample
        auto sample = new double[array_size]();
        random_sampling(data, sample, array_size);
        // calculate statistic from sample
        bin_stat[i] = operation(sample, array_size);
        delete[] sample;
    }
    // calculate the deviation of the statistic
    double result = stdev(bin_stat, bins);
    delete[] bin_stat;

    return result;
}

double mean(const double *array, int size) {
    double sum = std::accumulate(array, array + size, 0.0);
    return sum / size;
}

double stdev(const double *array, int size) {
    double result = 0;
    double mean_val = mean(array, size);
    for (int i = 0; i < size; i++) {
        result += pow(array[i] - mean_val, 2);
    }
    return sqrt(result / size);
}

int int_pow(int base, int exp) {
    int result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }
    return result;
}

double correlation(const double *data, int t, int t_max) {
    /// Calculates the correlation at time t
    /// This is done by multiplying the deviation of arr[0, t_max - t] with the deviation of arr[t, t_max]
    if (t > t_max)
        return 0;

    double result = 0;
    int sum_index = t_max - t;

    double mean_0 = mean(data, sum_index);
    double mean_t = mean(data + t, sum_index);

    //multiply the two deviations and return
    for (int t0 = 0; t0 < sum_index; t0++) {
        result += (data[t0] - mean_0) * (data[t0 + t] - mean_t);
    }
    return result / sum_index;
}

void autocorrelation(const double *data, double *result, int t_max) {
    /// Returns the autocorrelation function.
    /// Calculates the correlation at each t, from 0 up to t_max. At t=0, the two intervals are identical
    /// and the function is at a maximum; this is used as the normalisation such that the maximum of the
    /// autocorrelation curve is 1.

    //calculate the value for each t
    double norm = correlation(data, 0, t_max);
    for (int t = 0; t < t_max; t++) {
        result[t] = correlation(data, t, t_max) / norm;
    }
}
