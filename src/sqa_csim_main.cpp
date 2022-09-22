#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "sqa.hpp"
#include "sqa_csim_helper.hpp"

#ifndef U50
    #define U50 0
#endif
#ifndef WARP_MODE
    #define WARP_MODE 1
#endif

fp_t     Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM];
fpPack_t JcoupPack[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE];
fpPack_t JcoupPackBank0[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];
fpPack_t JcoupPackBank1[MAX_QUBIT_NUM][MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM];

u32_t       nTrotters = 16;
u32_t       nQubits   = MAX_QUBIT_NUM;
qubit_t     qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
qubitPack_t qubitsPack[MAX_TROTTER_NUM][MAX_QUBIT_NUM / PACKET_SIZE];
fp_t        h[MAX_QUBIT_NUM];
fp_t        prbNumbers[MAX_QUBIT_NUM];

const int  nSteps     = 100;     // default 500
const fp_t gammaStart = 10.0f;   // default 3.0f
const fp_t gammaEnd   = 0.001f;  // default 0.001f
const fp_t T          = 0.1f;    // default 0.3f

static void warp_execution(std::ofstream &energyLog);
static void unwarp_execution(std::ofstream &energyLog);

int main(int argc, char **argv)
{
    std::cout << "Current host settings:" << std::endl;
    std::cout << "* U50       : " << U50 << std::endl;
    std::cout << "* WARP_MODE : " << WARP_MODE << std::endl;
    std::cout << "* #SPIN     : " << MAX_QUBIT_NUM << std::endl;
    std::cout << "* #TROTTER  : " << nTrotters << std::endl;

    GenerateJcoupNP(nQubits, Jcoup, prbNumbers);
    GenerateHNP(nQubits, h);
    GenerateRandomQubitsState(nTrotters, nQubits, qubits);

#if U50
    PackQubits(qubitsPack, qubits);
    PackJcoup(JcoupPack, Jcoup);
    DispatchJcoup(JcoupPackBank0, JcoupPackBank1, JcoupPack);
#endif

    std::string   label = (U50) ? std::string("u50") : std::string("naive");
    std::ofstream energyLog(std::string("energy.") + label + std::string(".txt"));

#if WARP_MODE
    warp_execution(energyLog);
#else
    unwarp_execution(energyLog);
#endif
}

static void warp_execution(std::ofstream &energyLog)
{
    fp_t minEnergy = -1.0f;
    fp_t jperpMem[MAX_STEP_NUM];

    fp_t r        = pow(gammaEnd / gammaStart, 1.0f / (fp_t)(nSteps - 1));
    fp_t curGamma = gammaStart;
    for (u32_t i = 0; i < nSteps; i++) {
        fp_t jperp  = -0.5 * T * log(tanh(curGamma / (fp_t)nTrotters / T));
        jperpMem[i] = jperp;
        curGamma *= r;
    }

#if U50
    RunSQAHardware(nTrotters, nQubits, nSteps, 1.0f / T, 0, JcoupPackBank0, JcoupPackBank1, h,
                   jperpMem, qubitsPack);
#else
    RunSQASoftware(nTrotters, nQubits, nSteps, 1.0f / T, 0, Jcoup, h, jperpMem, qubits);
#endif

    for (u32_t i = 0; i < nSteps; i++) {
#if U50
        UnpackQubits(qubitsMemLogHW[i], qubits);
#endif
        for (u32_t t = 0; t < nTrotters; t++) {
#if U50
            fp_t energy = ComputeEnergyPerTrotter(nQubits, qubits[t], Jcoup, h);
#else
            fp_t energy = ComputeEnergyPerTrotter(nQubits, qubitsLogSW[i][t], Jcoup, h);
#endif
            energyLog << "T" << t << ": " << energy << std::endl;
            if (minEnergy < 0 || minEnergy > abs(energy)) minEnergy = abs(energy);
        }
        energyLog << "Energy: " << minEnergy << std::endl;
    }

    energyLog.close();
}

static void unwarp_execution(std::ofstream &energyLog)
{
    fp_t minEnergy = -1.0f;

    for (u32_t i = 0; i < nSteps; i++) {
        if ((i + 1) % 20 == 0) std::cout << (i + 1) << " iterations done..." << std::endl;

        fp_t gamma = gammaStart * (1.0f - ((fp_t)i / (fp_t)nSteps));
        fp_t jperp = -0.5 * T * log(tanh(gamma / (fp_t)nTrotters / T));

        fp_t rndNum[MAX_TROTTER_NUM][MAX_QUBIT_NUM];
        GenerateRandomNumber(nTrotters, rndNum, 1.0f / T);  // log(unif(rng)) / beta * 0.5f;

#if U50
        RunSQAHardwareOneStep(nTrotters, nQubits, jperp, h, rndNum, JcoupPackBank0, JcoupPackBank1,
                              qubitsPack);
        UnpackQubits(qubitsPack, qubits);
#else
        RunSQASoftwareOneStep(nTrotters, nQubits, jperp, h, rndNum, Jcoup, qubits);
#endif

        for (u32_t t = 0; t < nTrotters; t++) {
            fp_t energy = ComputeEnergyPerTrotter(nQubits, qubits[t], Jcoup, h);
            energyLog << "T" << t << ": " << energy << std::endl;
            if (minEnergy < 0 || minEnergy > abs(energy)) minEnergy = abs(energy);
        }

        energyLog << "Energy: " << minEnergy << std::endl;
    }

    energyLog.close();
}
