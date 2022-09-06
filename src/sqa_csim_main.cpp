#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "sqa.hpp"
#include "sqa_csim_helper.hpp"

#ifndef U50
    #define U50 1
#endif
#ifndef AM
    #define AM 0
#endif
#ifndef REPLAY
    #define REPLAY 0
#endif

// Jcoup
fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM];
fpPack_t Jcoup_pack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE];
fpPack_t Jcoup_pack_0[MAX_QUBIT_NUM]
                     [MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];
fpPack_t Jcoup_pack_1[MAX_QUBIT_NUM]
                     [MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];

int main(int argc, char **argv)
{
    // Dump value of macros
    std::cout << "Current host settings:" << std::endl;
    std::cout << "-> U50    : " << U50 << std::endl;
    std::cout << "-> AM     : " << AM << std::endl;
    std::cout << "-> REPLAY : " << REPLAY << std::endl;
    std::cout << "-> #SPIN  : " << MAX_QUBIT_NUM << std::endl;
    std::cout << "-> #TROT  : " << MAX_TROTTER_NUM << std::endl;

#if AM
    // Set nTrot and nSpin
    int nTrot = MAX_TROTTER_NUM;
    int nSpin = 18;
#else
    // Set nTrot and nSpin
    int nTrot = MAX_TROTTER_NUM;
    int nSpin = MAX_QUBIT_NUM;
#endif

    // Trotters
    qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
    qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE];

    // h
    fp_t h[MAX_QUBIT_NUM];

#if REPLAY
    // Read Trotters
    ReadRandomState(trotters, nTrot, nSpin, "../../../init_trotter.txt");
#else
    // Generate Trotters
    GenerateRandomState(trotters, nTrot, nSpin);
#endif

#if AM
    // Read Jcoup, h
    ReadJcoupAM(nSpin, Jcoup, "../../../J_r.txt");
    ReadHAM(nSpin, h, "../../../h_r.txt");
#else
    // Generate Jcoup, h
    fp_t rand_nums[MAX_QUBIT_NUM];
    GenerateJcoupNP(nSpin, Jcoup, rand_nums);
    GenerateHNP(nSpin, h);
#endif

#if U50
    // Convert into pack form
    PackTrotters(trotters_pack, trotters);
    PackJcoup(Jcoup_pack, Jcoup);
    // Dispatch
    for (u32_t i = 0; i < MAX_QUBIT_NUM; i++) {
        for (u32_t pack_ofst = 0;
             pack_ofst < MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM;
             pack_ofst++) {
            Jcoup_pack_0[i][pack_ofst] =
                Jcoup_pack[i][pack_ofst * JCOUP_BANK_NUM + 0];
            Jcoup_pack_1[i][pack_ofst] =
                Jcoup_pack[i][pack_ofst * JCOUP_BANK_NUM + 1];
        }
    }
#endif

    // Iteration Record
#if U50
    std::ofstream out_log("out_log.u50.txt");
    std::ofstream time_log("time_log.u50.txt");
#else
    std::ofstream out_log("out_log.basic.txt");
    std::ofstream time_log("time_log.basic.txt");
#endif

    // Iteration parameters
#if AM
    const int  iter        = 500;   // default 500
    const fp_t gamma_start = 3.0f;  // default 3.0f
    const fp_t T           = 0.3f;  // default 0.3f
#else
    const int     iter        = 500;     // default 500
    const fp_t    gamma_start = 3.0f;    // default 3.0f
    const fp_t    T           = 128.0f;  // default 0.3f
#endif

#if REPLAY
    // Read Log Random Number
    std::ifstream file_lrn("../../../log_rnd.txt");
#endif

    // Iteration of SQA
    for (int i = 0; i < iter; i++) {
        // Print progress
        if ((i + 1) % 20 == 0)
            std::cout << (i + 1) << " iterations done..." << std::endl;

        // Get jperp
        fp_t gamma = gamma_start * (1.0f - ((fp_t)i / (fp_t)iter));
        fp_t jperp = -0.5 * T * log(tanh(gamma / (fp_t)nTrot / T));

        // Generate Log Rand Number = log(unif(rng)) * nTrot
        fp_t log_rand_nums[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
#if REPLAY
        for (int t = 0; t < nTrot; t++) {
            for (int i = 0; i < nSpin; i++) { file_lrn >> log_rand_nums[t][i]; }
        }
#else
        GenerateLogRandomNumber(nTrot, log_rand_nums);
#endif

        // Set Timer
        std::chrono::system_clock::time_point start =
            std::chrono::system_clock::now();

#if U50
        // Run QMC-U50
        RunSQAHardware(trotters_pack, Jcoup_pack_0, Jcoup_pack_1, h, jperp,
                       1.0f / T, log_rand_nums);
        // Unpack Trotters
        UnpackTrotters(trotters_pack, trotters);
#else
        // Run QMC-Basic
        RunSQASoftware(nTrot, nSpin, trotters, Jcoup, h, jperp, 1.0f / T,
                                    log_rand_nums);
#endif

        // End Timer
        std::chrono::system_clock::time_point end =
            std::chrono::system_clock::now();
        time_log << "Run " << i << ": "
                 << std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                          start)
                        .count()
                 << " us" << std::endl;

        // Current state energy
        fp_t sum_energy = 0.0f;

        // Compute energy
        for (int t = 0; t < nTrot; t++) {
            // Get energy of each trotters
            fp_t energy = ComputeEnergyPerTrotter(nSpin, trotters[t], Jcoup, h);

            // Sum up
            sum_energy += energy;

            // Write out to log
            out_log << "T" << t << ": " << energy << std::endl;
        }

        // Write out total energy
        out_log << "Energy: " << sum_energy << std::endl;
    }

    // Close out_log
    out_log.close();
    time_log.close();

#if REPLAY
    // Close file_lrn
    file_lrn.close();
#endif
}
