#include <fstream>
#include <iostream>
#include <random>

#include "qmc_hy_driver.hpp"

#define PBSTR   "============================================================"
#define PBWIDTH 60

void PrintProgress(int run, int max_run)
{
    float percentage = (float)run / (float)max_run;
    int   val        = (int)(percentage * 100);
    int   lpad       = (int)(percentage * PBWIDTH);
    int   rpad       = PBWIDTH - lpad;
    printf("\r[STAT][%3d%%] [%.*s%*s] %4d/%4d", val, lpad, PBSTR, rpad, "", run,
           max_run);
    fflush(stdout);
}

int main(int argc, char *argv[])
{
    // Print usage
    if (argc != 6) {
        std::cout << "./host xclbin Jcoup h schedule n_trot" << std::endl;
        return -1;
    }

    bool  trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
    float jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM];
    float h[MAX_QUBIT_NUM];

    std::string jcoup_path = argv[2];
    std::string h_path     = argv[3];
    std::string sche_path  = argv[4];
    int         n_trot     = atoi(argv[5]);

    qmc_hy_driver qhd(argv[1], n_trot, jcoup, h);

    std::ifstream jcoup_stream(jcoup_path);
    std::ifstream h_stream(jcoup_path);
    std::ifstream sche_stream(sche_path);

    // Read jcoup and h
    for (int i = 0; i < MAX_QUBIT_NUM; i++) {
        h_stream >> h[i];
        for (int j = 0; j < MAX_QUBIT_NUM; j++) { jcoup_stream >> jcoup[i][j]; }
    }

    jcoup_stream.close();
    h_stream.close();

    // Iteration
    float         jperp, beta;
    std::ofstream out_stream("trotters.hy.txt");

    std::cout << "[INFO][!] -> Start Iterations of SQA" << std::endl;
    while (sche_stream >> jperp >> beta) {
        // Run
        qhd.run_qmc(jperp, beta);

        // Get Trotter
        qhd.get_trotters(trotters);

        // Write Trotter
        for (int t = 0; t < n_trot; t++) {
            for (int i = 0; i < MAX_QUBIT_NUM; i++) {
                out_stream << trotters[t][i] << " ";
            }
            out_stream << "\n";
        }
        out_stream << std::endl;
    }

    // Close files
    out_stream.close();
    sche_stream.close();
}