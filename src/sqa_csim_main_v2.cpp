#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "sqa.hpp"
#include "sqa_csim_helper.hpp"

#ifndef U50
    #define U50 1
#endif

fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM];
fpPack_t Jcoup_pack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE];
fpPack_t Jcoup_pack_0[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];
fpPack_t Jcoup_pack_1[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];

int main(int argc, char **argv)
{
    std::cout << "Current host settings:" << std::endl;
    std::cout << "* U50       : " << U50 << std::endl;
    std::cout << "* #SPIN     : " << MAX_QUBIT_NUM << std::endl;
    std::cout << "* #TROTTER  : " << MAX_TROTTER_NUM << std::endl;

    int         nTrot = MAX_TROTTER_NUM;
    int         nSpin = MAX_QUBIT_NUM;
    qubit_t     qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
    qubitPack_t qubitsPack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE];
    fp_t        h[MAX_QUBIT_NUM];
    fp_t        prbNumbers[MAX_QUBIT_NUM];

    GenerateJcoupNP(nSpin, Jcoup, prbNumbers);
    GenerateHNP(nSpin, h);
    GenerateRandomState(qubits, nTrot, nSpin);

#if U50
    PackTrotters(qubitsPack, qubits);
    PackJcoup(Jcoup_pack, Jcoup);
    for (u32_t i = 0; i < MAX_QUBIT_NUM; i++) {
        for (u32_t pack_ofst = 0; pack_ofst < MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM;
             pack_ofst++) {
            Jcoup_pack_0[i][pack_ofst] = Jcoup_pack[i][pack_ofst * JCOUP_BANK_NUM + 0];
            Jcoup_pack_1[i][pack_ofst] = Jcoup_pack[i][pack_ofst * JCOUP_BANK_NUM + 1];
        }
    }
#endif

#if U50
    std::ofstream out_log("out_log.u50.txt");
    std::ofstream time_log("time_log.u50.txt");
#else
    std::ofstream out_log("out_log.basic.txt");
    std::ofstream time_log("time_log.basic.txt");
#endif

    const int  iter        = 500;     // default 500
    const fp_t gamma_start = 3.0f;    // default 3.0f
    const fp_t T           = 128.0f;  // default 0.3f

    fp_t jperpMem[MAX_STEP_NUM];
    for (int i = 0; i < iter; i++) {
        fp_t gamma = gamma_start * (1.0f - ((fp_t)i / (fp_t)iter));
        fp_t jperp = -0.5 * T * log(tanh(gamma / (fp_t)nTrot / T));
        jperpMem[i] - jperp;
    }

#if U50
    RunSQAKernel(nTrot, nSpin, iter, 1.0f / T, 0, Jcoup_pack_0, Jcoup_pack_1, h, jperpMem,
                 qubitsPack);
#endif
}
