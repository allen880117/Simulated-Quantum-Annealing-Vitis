#pragma once
#include "sqa.hpp"

/* Progress Bar */
#define PBSTR   "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

/* Print the progress bar */
void PrintProgress(double percentage);

/* Read Random Initial State */
void ReadRandomState(qubit_t trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM], int nTrot, int nSpin,
                     std::string file_path);

/* Generate Random Initial State */
void GenerateRandomState(qubit_t trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM], int nTrot, int nSpin);

/* Generate Log Random Number = log(unif(rng)) * nTrot */
void GenerateLogRandomNumber(int nTrot, fp_t log_rand_nums[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                             fp_t beta, bool u50);

/* Compute Energy Summation of A Trotter*/
fp_t ComputeEnergyPerTrotter(int nSpin, qubit_t trotter[MAX_QUBIT_NUM],
                             fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t h[MAX_QUBIT_NUM]);

/* Convert trotters into pack form*/
void PackTrotters(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                  qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);

/* Unpack trotters */
void UnpackTrotters(qubitPack_t trotters_pack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
                    qubit_t     trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);

/* Convert Jcoup into pack form */
void PackJcoup(fpPack_t Jcoup_pack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE],
               fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM]);

/* Generate Jcoup of Number Partition */
void GenerateJcoupNP(int nSpin, fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                     fp_t rand_nums[MAX_QUBIT_NUM]);

/* Read Jcoup of AM */
void ReadJcoupAM(int nSpin, fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], std::string J_file_path);

/* Generate H of Number Partition */
void GenerateHNP(int nSpin, fp_t h[MAX_QUBIT_NUM]);

/* Read H of AM */
void ReadHAM(int nSpin, fp_t h[MAX_QUBIT_NUM], std::string h_file_path);
