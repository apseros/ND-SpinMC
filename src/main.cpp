#include "metropolis.h"
#include "arguments.h"

int main(int argc, char *argv[]) {
    Arguments settings = Arguments(argc, argv);
    int m = settings.get_model();
    if (m==0) {Ising I = Ising(settings);}
    else if (m==1) {XY_model I = XY_model(settings);}
    else if (m==2) {Heisenberg I = Heisenberg(settings);}
    else {std::cout << "Please specify model (0, 1, 2)" << std::endl; return 0;}
    return 0;
}