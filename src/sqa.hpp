#pragma once
#include "ap_int.h"
#include "hls_math.h"

#define PRAGMA_SUB(PRAG) _Pragma(#PRAG)
#define CTX_PRAGMA(PRAG) PRAGMA_SUB(PRAG)

#if __SYNTHESIS__
    #define USING_STD_RNG 0
#else
    #define USING_STD_RNG 0
#endif

#if USING_STD_RNG
    #include <random>
#endif

#define RESOURCE_CONSTRAINT 1

/* Common parameters */
#if __SYNTHESIS__
    #define MAX_TROTTER_NUM     16
    #define MAX_QUBIT_NUM       4096
    #define MAX_STEP_NUM        1024
    #define PACKET_SIZE         64
    #define LOG2_PACKET_SIZE    6
    #define JCOUP_BANK_NUM      2
    #define LOG2_JCOUP_BANK_NUM 1
    #define NUM_FADD            64
#else
    #define MAX_TROTTER_NUM     4
    #define MAX_QUBIT_NUM       32
    #define MAX_STEP_NUM        1024
    #define PACKET_SIZE         16
    #define LOG2_PACKET_SIZE    4
    #define JCOUP_BANK_NUM      2
    #define LOG2_JCOUP_BANK_NUM 1
    #define NUM_FADD            64
#endif

/* Generic Type */
typedef unsigned int u32_t;
typedef int          i32_t;
typedef float        fp_t;
typedef ap_uint<1>   qubit_t;

/* Aggreate Type */
typedef struct {
    fp_t data[PACKET_SIZE];
} fpPack_t;
typedef ap_uint<PACKET_SIZE> qubitPack_t;

/* Column number of each cache or memory bank */
#define NUM_COL_QUBIT_MEM      (MAX_QUBIT_NUM / PACKET_SIZE)
#define NUM_COL_QUBIT_CACHE    (NUM_COL_QUBIT_MEM)
#define NUM_COL_JCOUP_MEM_BANK (MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM)
#define NUM_COL_RNDNUM_CACHE   (MAX_QUBIT_NUM)

/* Prototype of Functions */
extern "C" {
void RunSQAHardwareOneStep(u32_t nTrotters, u32_t nQubits, fp_t jperp, fp_t hCache[MAX_QUBIT_NUM],
                           fp_t        rndNumMem[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE],
                           fpPack_t    jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                           fpPack_t    jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                           qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE]);
void RunSQAHardware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                    fpPack_t jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fpPack_t jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fp_t hMem[MAX_QUBIT_NUM], fp_t jperpMem[MAX_STEP_NUM],
                    qubitPack_t qubitsMem[MAX_TROTTER_NUM][NUM_COL_QUBIT_MEM],
                    qubitPack_t qubitsHistory[MAX_STEP_NUM][MAX_TROTTER_NUM][NUM_COL_QUBIT_MEM]);
}

void RunSQASoftwareOneStep(u32_t nTrotters, u32_t nQubits, fp_t jperp, fp_t h[MAX_QUBIT_NUM],
                           fp_t    rndNum[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                           fp_t    jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                           qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);
void RunSQASoftware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                    fp_t jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t hMem[MAX_QUBIT_NUM],
                    fp_t jperpMem[MAX_STEP_NUM], qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                    qubit_t qubitsHistory[MAX_STEP_NUM][MAX_TROTTER_NUM][MAX_QUBIT_NUM]);
