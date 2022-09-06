#pragma once
#include "hls_math.h"
#include "hls_stream.h"

#define PRAGMA_SUB(PRAG) _Pragma(#PRAG)
#define CTX_PRAGMA(PRAG) PRAGMA_SUB(PRAG)

/* Common parameters */
#if __SYNTHESIS__
    #define MAX_TROTTER_NUM     4
    #define MAX_QUBIT_NUM       4096
    #define PACKET_SIZE         64
    #define LOG2_PACKET_SIZE    6
    #define JCOUP_BANK_NUM      2
    #define LOG2_JCOUP_BANK_NUM 1
    #define NUM_FADD            64
#else
    #define MAX_TROTTER_NUM     4
    #define MAX_QUBIT_NUM       32
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
#define NUM_COL_QUBIT_CACHE    (MAX_QUBIT_NUM / PACKET_SIZE)
#define NUM_COL_JCOUP_MEM_BANK (MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM)
#define NUM_COL_RNDNUM_MEM     (MAX_QUBIT_NUM)

/* Prototype of Functions */
extern "C" {
void RunSQAHardware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t jperp,
                    fp_t hCache[MAX_QUBIT_NUM], fp_t rndNumMem[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM],
                    fpPack_t    jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fpPack_t    jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE]);
}

void RunSQASoftware(const int nTrot, const int nSpin,
                    qubit_t    trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                    const fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], const fp_t h[MAX_QUBIT_NUM],
                    const fp_t jperp, const fp_t Beta,
                    const fp_t logRandNumber[MAX_TROTTER_NUM][MAX_QUBIT_NUM]);
