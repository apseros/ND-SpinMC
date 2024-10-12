#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <iostream>
#include <unistd.h>
#include <libgen.h>
#include "metropolis.h"

#define OPT_ARGS "L:d:N:m:t:B:x:A:w:p:h"
#define DEFAULT_LENGTH        32
#define DEFAULT_DIMENSION     2
#define DEFAULT_ITERATIONS    5000
#define DEFAULT_B_FIELD		  0
#define DEFAULT_WOLFF		  0
#define DEFAULT_EXCHANGE      1
#define DEFAULT_ANISOTROPY    0.01
#define DEFAULT_TEMPERATURE   0
#define DEFAULT_MODEL         1
#define DEFAULT_PRINTS        10

class Arguments {
protected:
    std::string prog;
    int length = DEFAULT_LENGTH, dim = DEFAULT_DIMENSION, iterations = DEFAULT_ITERATIONS, model = DEFAULT_MODEL,
            intermediate_print = DEFAULT_PRINTS;
    double exchange = DEFAULT_EXCHANGE, anisotropy = DEFAULT_ANISOTROPY, magnetic_field = DEFAULT_B_FIELD,
            temperature = DEFAULT_TEMPERATURE;
    bool wolff;
    void validate();
    friend class Metropolis;
public:
    Arguments() = delete;
    Arguments(int argc, char* argv[]);
    void usage();
    [[nodiscard]] int get_model() const;
    void print() const;
};

#endif
