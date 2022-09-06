#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "sqa.hpp"
#include "sqa_csim_helper.hpp"

#ifndef U50
    #define U50 0
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

    fp_t minEnergy = -1.0f;

    for (int i = 0; i < iter; i++) {
        if ((i + 1) % 20 == 0) std::cout << (i + 1) << " iterations done..." << std::endl;
        fp_t gamma = gamma_start * (1.0f - ((fp_t)i / (fp_t)iter));
        fp_t jperp = -0.5 * T * log(tanh(gamma / (fp_t)nTrot / T));
        fp_t log_rand_nums[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
        GenerateLogRandomNumber(nTrot, log_rand_nums, 1.0f / T, (U50 == 1));

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
#if U50
        RunSQAHardware(nTrot, nSpin, jperp, h, log_rand_nums, Jcoup_pack_0, Jcoup_pack_1,
                       qubitsPack);
        UnpackTrotters(qubitsPack, qubits);
#else
        RunSQASoftware(nTrot, nSpin, qubits, Jcoup, h, jperp, 1.0f / T, log_rand_nums);
#endif
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        time_log << "Run " << i << ": "
                 << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
                 << " us" << std::endl;

        for (int t = 0; t < nTrot; t++) {
            fp_t energy = ComputeEnergyPerTrotter(nSpin, qubits[t], Jcoup, h);
            out_log << "T" << t << ": " << energy << std::endl;
            if (minEnergy < 0 || minEnergy > abs(energy)) minEnergy = abs(energy);
        }
        out_log << "Energy: " << minEnergy << std::endl;
    }

    // Close out_log
    out_log.close();
    time_log.close();
}
