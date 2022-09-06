#ifndef _SQA_HPP_
#define _SQA_HPP_

/* _QMC_U50_R4_HPP */

#include "hls_math.h"
#include "hls_stream.h"

#define PRAGMA_SUB(PRAG) _Pragma(#PRAG)
#define CTX_PRAGMA(PRAG) PRAGMA_SUB(PRAG)

#if __SYNTHESIS__
    #define NUM_TROT         4
    #define NUM_SPIN         4096
    #define PACKET_SIZE      64
    #define LOG2_PACKET_SIZE 6
    #define NUM_STREAM       2
    #define LOG2_NUM_STREAM  1
    #define NUM_FADD         64
#else
    #define NUM_TROT         4
    #define NUM_SPIN         32
    #define PACKET_SIZE      16
    #define LOG2_PACKET_SIZE 4
    #define NUM_STREAM       2
    #define LOG2_NUM_STREAM  1
    #define NUM_FADD         64
#endif

typedef unsigned int u32_t;
typedef int          i32_t;

typedef float fp_t;
typedef struct {
    fp_t data[PACKET_SIZE];
} fp_pack_t;

typedef ap_uint<1>           spin_t;
typedef ap_uint<PACKET_SIZE> spin_pack_t;

extern "C" {
void RunSQAHardware(spin_pack_t     trotters[NUM_TROT][NUM_SPIN / PACKET_SIZE],
                    const fp_pack_t Jcoup_0[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
                    const fp_pack_t Jcoup_1[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
                    const fp_t h[NUM_SPIN], const fp_t Jperp, const fp_t Beta,
                    const fp_t logRand[NUM_TROT][NUM_SPIN]);
}

void RunSQASoftware(const int nTrot, const int nSpin, spin_t trotters[NUM_TROT][NUM_SPIN],
                    const fp_t Jcoup[NUM_SPIN][NUM_SPIN], const fp_t h[NUM_SPIN], const fp_t Jperp,
                    const fp_t Beta, const fp_t logRandNumber[NUM_TROT][NUM_SPIN]);

#endif