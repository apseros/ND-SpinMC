#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <mutex>
#include "arguments.h"
void init_nn_map_helical(int*, int, int);

class Metropolis {
protected:
    int length{}, iterations{}, dimensions{}, volume{}, exchange{}, print_steps{};
    int *nn_map{};
    double magnetic_field{}, temperature{}, anisotropy{};
    bool wolff{};
    void metropolis_step();
    virtual void propose_change(int) {};
    virtual void initialise_spins() {};
    virtual void wolff_step();
    virtual void calculate_order(int) {};
public:
    Metropolis() = delete;
    Metropolis(int n_len, int n_iter, int n_dim, int ex, double field, double temp, bool);
    explicit Metropolis(class Arguments& a);
    ~Metropolis();
    void run();
    virtual void print_spins() {};
};

class Ising : public Metropolis{
protected:
    int *spins, *magnetisation;
    void propose_change(int) override;
    void initialise_spins() override;
    void calculate_order(int) override;
    void wolff_step() override;
public:
    [[maybe_unused]] Ising(int , int , int , int , double , double , bool);
    explicit Ising(class Arguments &a);
    ~Ising();
    void print_spins() override;
};

#define ROTATION_FRACTION 0.05
class XY_model : public Metropolis{
protected:
    double *spins, *magnetisation;
    void propose_change(int) override;
    void initialise_spins() override;
    void calculate_order(int) override;
    void wolff_step() override;
public:
    [[maybe_unused]] XY_model(int , int , int , int , double , double , bool);
    explicit XY_model(class Arguments &a);
    ~XY_model();
    void print_spins() override;
};

class Heisenberg : public Metropolis{
protected:
    double *spin_theta, *spin_phi, *magnetisation;
    void propose_change(int) override;
    void initialise_spins() override;
    void calculate_order(int) override;
    void wolff_step() override;
public:
    [[maybe_unused]] [[maybe_unused]] Heisenberg(int , int , int , int , double , double , bool);
    explicit Heisenberg(class Arguments &a);
    ~Heisenberg();
    void print_spins() override;
};

#endif
