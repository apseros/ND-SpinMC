#include "metropolis.h"
#include "utilities.h"
#include "arguments.h"
#include <iostream>
#include <thread>
#include <cmath>
#include <iomanip>

Metropolis::Metropolis(int n_len, int n_iter, int n_dim, int ex, double field, double temp, bool wolff)
        : length(n_len), iterations(n_iter), dimensions(n_dim), exchange(ex), magnetic_field(field), temperature(temp),
          wolff(wolff) {
    volume = int_pow(length, dimensions);
    nn_map = new int[dimensions * 2 * volume]();
    init_nn_map_helical(nn_map, length, dimensions);
}

Metropolis::Metropolis(class Arguments &a) : length(a.length), iterations(a.iterations), dimensions(a.dim),
                                             exchange(static_cast<int>(a.exchange)), magnetic_field(a.magnetic_field),
                                             temperature(a.temperature), wolff(a.wolff), anisotropy(a.anisotropy) {
    volume = int_pow(length, dimensions);
    nn_map = new int[dimensions * 2 * volume]();
    print_steps = iterations / a.intermediate_print;
    init_nn_map_helical(nn_map, length, dimensions);
}

Metropolis::~Metropolis() {
    delete[] nn_map;
}

void Metropolis::run() {
    int n_prints = print_steps;
    for (int i = 0; i < iterations; i++) {
        wolff ? wolff_step() : metropolis_step();
        calculate_order(i);
        if (i == n_prints - 1) {
            print_spins();
            if (print_steps < iterations) {
                std::cout << std::endl;
                n_prints += print_steps;
            }
        }
    }
}


void Metropolis::metropolis_step() {
    std::uniform_int_distribution<int> dist_idx{0, volume - 1};
    for (int j = 0; j < volume; j++) {
        int rand_idx = dist_idx(rng);
        propose_change(rand_idx);
    }
}

void Metropolis::wolff_step() {
    wolff = false;
    metropolis_step();
}


void init_nn_map_helical(int *nn_map, int length, int dim) {
    /// helical boundary conditions imply helicity in n-1 dimensions and cyclicity in 1 dimension
    /// helicity -> cyclicity as N -> infinity
    /// since in N dimensions, each lattice point has 2N nearest neighbours, these nearest neighbours are explored
    /// in N dimensions, the nearest neighbours have the array indices (relative to the current point):
    /// +1, -1; +L, -L; +L^2, -L^2 ; ... ; -L^(N-1), +L^(N-1).
    int vol = int_pow(length, dim);
    int temp, temp_dim_power;
    for (int i = 0; i < vol; i++) {
        for (int j = 0; j < dim; j++) {
            temp_dim_power = int_pow(length, j);
            //positive direction
            temp = i + temp_dim_power;
            nn_map[i * 2 * dim + 2 * j] = (temp > vol ? temp - vol : temp);
            //negative direction
            temp = i - temp_dim_power;
            nn_map[i * 2 * dim + 2 * j + 1] = (temp < 0 ? temp + vol : temp);
        }
    }
}

void Ising::initialise_spins() {
    std::uniform_int_distribution<int> distrib{0, 1};
    for (int i = 0; i < volume; i++) {
        spins[i] = (2 * distrib(rng)) - 1;
    }
}

[[maybe_unused]] Ising::Ising(int n_len, int n_iter, int n_dim, int ex, double field, double temp, bool wolff)
        : Metropolis(n_len, n_iter, n_dim, ex, field, temp, wolff) {
    spins = new int[volume];
    Ising::initialise_spins();
    magnetisation = new int[iterations]();
    run();
}

Ising::Ising(Arguments &a) : Metropolis(a) {
    spins = new int[volume];
    Ising::initialise_spins();
    magnetisation = new int[iterations]();
    run();
}


void Ising::print_spins() {
    for (int i = 0; i < volume; i++) {
        if (i % length == 0 && i != 0)
            std::cout << std::endl;
        if (spins[i] == -1)
            std::cout << ".";
        else
            std::cout << "0";
    }
    std::cout << std::endl;
}


void Ising::propose_change(int index) {
    std::uniform_real_distribution<double> fluctuation{0.0, 1.0};
    double delta_energy = 0;

    for (int i = 0; i < 2 * dimensions; i++) {
        delta_energy += spins[nn_map[index * 2 * dimensions + i]];
    }
    delta_energy *= 2 * exchange * spins[index];
    delta_energy += magnetic_field * spins[index];

    if (delta_energy <= 0 or fluctuation(rng) < exp(-delta_energy / temperature)) spins[index] *= -1;
}

Ising::~Ising() {
    delete[] spins;
    delete[] magnetisation;
}

void Ising::calculate_order(int iteration) {
    magnetisation[iteration] = std::accumulate(spins, spins + volume, 0);
}


void XY_model::initialise_spins() {
    std::uniform_real_distribution<double> distrib{0, 2 * M_PI};
    for (int i = 0; i < volume; i++) spins[i] = distrib(rng);
}

[[maybe_unused]] XY_model::XY_model(int n_len, int n_iter, int n_dim, int ex, double field, double temp, bool wolff)
        : Metropolis(n_len, n_iter, n_dim, ex, field, temp, wolff) {
    spins = new double[volume];
    XY_model::initialise_spins();
    magnetisation = new double[iterations * 2]();
    run();
}

XY_model::XY_model(Arguments &a) : Metropolis(a) {
    spins = new double[volume];
    XY_model::initialise_spins();
    magnetisation = new double[iterations * 2]();
    run();
}

void XY_model::print_spins() {
    for (int i = 0; i < volume; i++) {
        if (i % length == 0 && i != 0)
            std::cout << std::endl;
        std::cout << std::fixed;
        std::cout << std::setprecision(2);
        std::cout << spins[i] << " ";
    }
    std::cout << std::endl;
}

void XY_model::propose_change(int index) {
    std::uniform_real_distribution<double> delta_phi_dist{-ROTATION_FRACTION, ROTATION_FRACTION};
    std::uniform_real_distribution<double> fluctuation{0.0, 1.0};

    double delta_energy = 0;
    double central_spin = spins[index];
    double delta_phi = delta_phi_dist(rng);

    for (int i = 0; i < 2 * dimensions; i++) {
        double current_site = spins[nn_map[index * 2 * dimensions + i]];
        delta_energy -=
                exchange * 2 * (cos(central_spin + delta_phi - current_site) - cos(central_spin - current_site));
        delta_energy -= anisotropy * (cos(3 * (central_spin + delta_phi)) - cos(3 * central_spin));
    }

    if (delta_energy <= 0 or fluctuation(rng) < exp(-delta_energy / temperature)) spins[index] += delta_phi;
}

XY_model::~XY_model() {
    delete[] spins;
    delete[] magnetisation;
}


void XY_model::calculate_order(int iteration) {
    double x_comp = 0;
    double y_comp = 0;
    for (int i = 0; i < volume; i++) {
        x_comp += cos(spins[i]);
        y_comp += sin(spins[i]);
    }
    magnetisation[iteration * 2] = x_comp;
    magnetisation[iteration * 2 + 1] = y_comp;
}

void XY_model::wolff_step() {
    Metropolis::wolff_step();
}


void Heisenberg::initialise_spins() {
    std::uniform_real_distribution<double> distrib{0, 2 * M_PI};
    for (int i = 0; i < volume; i++) {
        spin_theta[i] = distrib(rng);
        spin_phi[i] = distrib(rng);
    }
}

[[maybe_unused]] Heisenberg::Heisenberg(int n_len, int n_iter, int n_dim, int ex, double field, double temp, bool wolff)
        : Metropolis(n_len, n_iter, n_dim, ex, field, temp, wolff) {
    spin_theta = new double[volume];
    spin_phi = new double[volume];
    Heisenberg::initialise_spins();

    magnetisation = new double[iterations * 3]();
    run();
}

Heisenberg::Heisenberg(Arguments &a) : Metropolis(a) {
    spin_theta = new double[volume];
    spin_phi = new double[volume];
    Heisenberg::initialise_spins();

    magnetisation = new double[iterations * 3]();
    run();
}


void Heisenberg::propose_change(int index) {
    std::uniform_real_distribution<double> delta_phi_dist{-ROTATION_FRACTION, ROTATION_FRACTION};
    std::uniform_real_distribution<double> fluctuation{0.0, 1.0};

    double delta_energy = 0;
    double t = spin_theta[index];
    double p = spin_phi[index];

    double dp = delta_phi_dist(rng);
    double dt = delta_phi_dist(rng);

    for (int i = 0; i < 2 * dimensions; i++) {
        double tt = spin_theta[nn_map[index * 2 * dimensions + i]];
        double pp = spin_phi[nn_map[index * 2 * dimensions + i]];
        delta_energy -= exchange * (sin(t + dt) * cos(p + dp) - sin(t) * cos(p)) * sin(tt) * cos(pp);
        delta_energy -= exchange * (sin(t + dt) * sin(p + dp) - sin(t) * sin(p)) * sin(tt) * sin(pp);
        delta_energy -= exchange * (cos(t + dt) - cos(t)) * cos(tt);
    }

    if (delta_energy <= 0 or fluctuation(rng) < exp(-delta_energy / temperature)) {
        spin_theta[index] += dt;
        spin_phi[index] += dp;
    }
    spin_theta[index] = fmod(spin_theta[index], 2 * M_PI);
    spin_phi[index] = fmod(spin_phi[index], 2 * M_PI);
}

Heisenberg::~Heisenberg() {
    delete[] spin_theta;
    delete[] spin_phi;
    delete[] magnetisation;
}


void Heisenberg::calculate_order(int iteration) {
    double x_comp = 0;
    double y_comp = 0;
    double z_comp = 0;
    for (int i = 0; i < volume; i++) {
        x_comp += sin(spin_theta[i]) * cos(spin_phi[i]);
        y_comp += sin(spin_theta[i]) * sin(spin_phi[i]);
        z_comp += cos(spin_theta[i]);
    }
    magnetisation[iteration * 3] = x_comp;
    magnetisation[iteration * 3 + 1] = y_comp;
    magnetisation[iteration * 3 + 2] = z_comp;
}


void Heisenberg::print_spins() {
    for (int i = 0; i < volume; i++) {
        std::cout << spin_theta[i] << " " << spin_phi[i] << std::endl;
    }
}

void Heisenberg::wolff_step() {
    Metropolis::wolff_step();
}


void Ising::wolff_step() {
    //this function executes a single step of the wolff algorithm, changing the
    //S[] array used in the arguments and return the new magnetisation using the
    //mag argument.

    //define the stack of lattice points in the cluster
    std::vector<int> stack;
    std::uniform_real_distribution<double> fluctuation{0.0, 1.0};
    std::uniform_int_distribution<int> dist_idx{0, volume - 1};

    int index, cluster_size;
    auto flag_array = new bool[volume]();

    //randomly choose a seed, flag it, add it to the cluster and increase cluster size
    int cluster_seed = dist_idx(rng);
    flag_array[cluster_seed] = true;
    stack.push_back(cluster_seed);
    cluster_size = 1;

    //define the probability of adding a new element to the cluster and the pseudo-random
    //number which will be used for randomness, x
    double prob = 1 - exp(-2 / temperature); // move to delta energy

    //while the stack of lattice points of which their surrounding will be explored evaluate:
    while (!stack.empty()) {     //repeat until there are no more points in the stack
        //get the last lattice point from the stack and remove it from the stack
        index = stack.back();
        stack.pop_back();

        //these are evaluated using a for loop
        for (int i = 0; i < 2 * dimensions; i++) {
            //get a nearest neighbour from the NNmap array
            int temp = nn_map[index * 2 * dimensions + i];
            if (!flag_array[temp]) {
                //randomly add the new point if it is also pointing in the same direction
                //as the seed spins
                if (fluctuation(rng) < prob && spins[temp] == exchange * spins[index]) {
                    //if these tests are successfully, flag the new lattice point, add it to the stack
                    //and increase the lattice size
                    flag_array[temp] = true;
                    stack.push_back(temp);
                    cluster_size++;
                }
            }
        }
    }

    //flip all members of the cluster
    for (int k = 0; k < volume; k++) {
        if (flag_array[k]) {
            spins[k] *= -1;
        }
    }

    //delete dynamic arrays and return the magnetisation
    delete[] flag_array;
}