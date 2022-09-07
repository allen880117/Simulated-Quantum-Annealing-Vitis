#include "sqa_csim_helper.hpp"

#include <cmath>
#include <fstream>
#include <random>
#include <string>

void PrintProgress(double percentage)
{
    int val  = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
    if (val == 100) printf("\n");
}

void GenerateRandomQubitsState(int nTrot, int nSpin,
                               qubit_t trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    /* Random Device & Random Number Generator */
    static std::mt19937                       rng(12345);
    static std::uniform_int_distribution<int> unif(0, 1);
    /* Generate */
    for (int i = 0; i < nTrot; i++) {
        for (int j = 0; j < nSpin; j++) { trotters[i][j] = (bool)unif(rng); }
    }
}

void GenerateRandomNumber(int nTrot, fp_t log_rand_nums[MAX_TROTTER_NUM][MAX_QUBIT_NUM], fp_t beta)
{
    // Random Generator
    static std::mt19937                  rng(12347);
    std::uniform_real_distribution<fp_t> unif(0.0, 1.0);

    // Generate Log Rand Number = log(unif(rng)) * nTrot
    for (int t = 0; t < MAX_TROTTER_NUM; t++) {
        for (int i = 0; i < MAX_QUBIT_NUM; i++) {
            log_rand_nums[t][i] = log(unif(rng)) / beta * 0.5f;
        }
    }
}

fp_t ComputeEnergyPerTrotter(int nSpin, qubit_t trotter[MAX_QUBIT_NUM],
                             fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t h[MAX_QUBIT_NUM])
{
    /* Follow the Formula */

    /* Total Energy */
    fp_t H = 0;

    /* Traverse All Spins */
    for (int i = 0; i < nSpin; i++) {
        /* Two-Spin Energy */
        for (int j = 0; j < nSpin; j++) {
            if (trotter[i] == trotter[j])
                H += Jcoup[i][j];
            else
                H -= Jcoup[i][j];
        }
        /* One-Spin Energy */
        if (trotter[i])
            H += h[i];
        else
            H -= h[i];
    }

    /* Return */
    return (H);
}

void PackQubits(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    for (int t = 0; t < MAX_TROTTER_NUM; t++) {
        for (int i = 0; i < MAX_QUBIT_NUM / PACKET_SIZE; i++) {
            for (int k = 0; k < PACKET_SIZE; k++) {
                trotters_pack[t][i][k] = trotters[t][i * PACKET_SIZE + k];
            }
        }
    }
}

void UnpackQubits(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                  qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    for (int t = 0; t < MAX_TROTTER_NUM; t++) {
        for (int i = 0; i < MAX_QUBIT_NUM / PACKET_SIZE; i++) {
            for (int k = 0; k < PACKET_SIZE; k++) {
                trotters[t][i * PACKET_SIZE + k] = trotters_pack[t][i][k];
            }
        }
    }
}

void PackJcoup(fpPack_t JcoupPack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
               fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM])
{
    for (int i = 0; i < MAX_QUBIT_NUM; i++) {
        for (int j = 0; j < MAX_QUBIT_NUM / PACKET_SIZE; j++) {
            for (int k = 0; k < PACKET_SIZE; k++) {
                JcoupPack[i][j].data[k] = Jcoup[i][j * PACKET_SIZE + k];
            }
        }
    }
}

void DispatchJcoup(
    fpPack_t JcoupPackBank0[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM],
    fpPack_t JcoupPackBank1[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM],
    fpPack_t JcoupPack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE])
{
    for (u32_t i = 0; i < MAX_QUBIT_NUM; i++) {
        for (u32_t pack_ofst = 0; pack_ofst < MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM;
             pack_ofst++) {
            JcoupPackBank0[i][pack_ofst] = JcoupPack[i][pack_ofst * JCOUP_BANK_NUM + 0];
            JcoupPackBank1[i][pack_ofst] = JcoupPack[i][pack_ofst * JCOUP_BANK_NUM + 1];
        }
    }
}

void GenerateJcoupNP(int nSpin, fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                     fp_t rand_nums[MAX_QUBIT_NUM])
{
    for (int i = 0; i < nSpin; i++) { rand_nums[i] = (float)(i + 1) / (float)nSpin; }

    for (int i = 0; i < nSpin; i++) {
        for (int j = 0; j < nSpin; j++) { Jcoup[i][j] = rand_nums[i] * rand_nums[j]; }
    }
}

void GenerateHNP(int nSpin, fp_t h[MAX_QUBIT_NUM])
{
    for (int i = 0; i < nSpin; i++) { h[i] = 0.0f; }
}
