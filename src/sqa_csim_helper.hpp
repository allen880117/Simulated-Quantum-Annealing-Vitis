#pragma once
#include "sqa.hpp"

/* Progress Bar */
#define PBSTR   "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

/* Print the progress bar */
void PrintProgress(double percentage);

/* Generate Random Initial State */
void GenerateRandomQubitsState(int nTrot, int nSpin,
                               qubit_t trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);

/* Generate Log Random Number = log(unif(rng)) * nTrot */
void GenerateRandomNumber(int nTrot, fp_t log_rand_nums[MAX_TROTTER_NUM][MAX_QUBIT_NUM], fp_t beta);

/* Compute Energy Summation of A Trotter*/
fp_t ComputeEnergyPerTrotter(int nSpin, qubit_t trotter[MAX_QUBIT_NUM],
                             fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t h[MAX_QUBIT_NUM]);

/* Convert trotters into pack form*/
void PackQubits(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);

/* Unpack trotters */
void UnpackQubits(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                  qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);

/* Convert Jcoup into pack form */
void PackJcoup(fpPack_t JcoupPack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
               fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM]);

/* Dispatch Jcoup into multile banks */
void DispatchJcoup(
    fpPack_t JcoupPackBank0[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM],
    fpPack_t JcoupPackBank1[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM],
    fpPack_t JcoupPack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE]);

/* Generate Jcoup of Number Partition */
void GenerateJcoupNP(int nSpin, fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                     fp_t rand_nums[MAX_QUBIT_NUM]);

/* Generate H of Number Partition */
void GenerateHNP(int nSpin, fp_t h[MAX_QUBIT_NUM]);
