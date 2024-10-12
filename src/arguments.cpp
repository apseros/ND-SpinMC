#include "arguments.h"

Arguments::Arguments(int argc, char **argv) {
    // function to obtain command line arguments and assign them to class variables
    // if something is unrecognised, this function calls the usage function
    int c;
    bool error_flg = false;
    prog.assign((char *) basename(argv[0]));
    while (!error_flg && (c = getopt(argc, argv, OPT_ARGS)) != -1) {
        switch (c) {
            case 'L':
                length = std::stoi(optarg);
                break;
            case 'd':
                dim = std::stoi(optarg);
                break;
            case 'N':
                iterations = std::stoi(optarg);
                break;
            case 't':
                temperature = std::stod(optarg);
                break;
            case 'm':
                model = std::stoi(optarg);
                break;
            case 'B':
                magnetic_field = std::stod(optarg);
                break;
            case 'w':
                wolff = std::stoi(optarg);
                break;
            case 'x':
                exchange = std::stod(optarg);
                break;
            case 'A':
                anisotropy = std::stod(optarg);
                break;
            case 'p':
                intermediate_print = std::stoi(optarg);
                break;
            case 'h':
                usage();
                break;
            default:
                std::cout << "Unknown option: " << static_cast<char>(c) << std::endl;
                error_flg = true;
        }
        if (error_flg) usage();
    }
    validate();
}

void Arguments::usage() {
    //prints the command line options for this program and their default
    //values and then exits with a warning.
    std::cout << "Usage: " << prog << "[options]" << std::endl
              << "-L: Lattice length (default = " << DEFAULT_LENGTH << ")" << std::endl
              << "-d: Lattice dimension (default = " << DEFAULT_DIMENSION << ")" << std::endl
              << "-N: Number of thermalization iterations (default = " << DEFAULT_ITERATIONS << ")" << std::endl
              << "-t: Temperature (default = " << DEFAULT_TEMPERATURE << ")" << std::endl
              << "-m: Model 0 for Ising, 1 for XY, 2 for Heisenberg (default = " << DEFAULT_MODEL << ")" << std::endl
              << "-B: Magnetic Field strength (default = " << DEFAULT_B_FIELD << ")" << std::endl
              << "-w: Wolff method, 0 for Metropolis, 1 for Wolff clusters (default = " << DEFAULT_WOLFF << ")"
              << std::endl
              << "-A: Anisotropy: (default = " << DEFAULT_ANISOTROPY << ")" << std::endl
              << "-x: Exchange interaction, (default = " << DEFAULT_EXCHANGE << ")" << std::endl
              << "-p: Number of intermediate prints of the configuration, (default = " << DEFAULT_PRINTS << ")"
              << std::endl
              << "Monte-Carlo simulation of the N-dimensional Ising/XY/Heisenberg models using" << std::endl
              << "the Metropolis algorithm. Wolff clustering is also an option for the Ising model." << std::endl
              << "The program can be used to print a configuration after -N thermalization steps" << std::endl;
    exit(1);
}

void Arguments::validate() {
    if (dim < 1) dim = DEFAULT_DIMENSION;
    if (length < 2) length = DEFAULT_LENGTH;
    if (iterations < 1) iterations = DEFAULT_ITERATIONS;
    if (temperature < 0) temperature = DEFAULT_TEMPERATURE;
    if (model < 0 || model > 2) model = DEFAULT_MODEL;
    if (intermediate_print < 0) intermediate_print = DEFAULT_PRINTS;
}

void Arguments::print() const {
    std::cout << "Current Parameters: " << std::endl
              << "Lattice length: " << length << std::endl
              << "Lattice dimension: " << dim << std::endl
              << "Number of thermalization iterations: " << iterations << std::endl
              << "Model: " << model << std::endl
              << "Temperature: " << temperature << std::endl
              << "Magnetic Field strength: " << magnetic_field << std::endl
              << "Wolff method: " << wolff << std::endl
              << "Anisotropy: " << anisotropy << std::endl
              << "Intermediate prints: " << intermediate_print << std::endl
              << "Exchange interaction: " << exchange << std::endl;
}


int Arguments::get_model() const {
    return model;
}