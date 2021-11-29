#ifndef _SQA_HPP_
#define _SQA_HPP_

#include "hls_math.h"
#include "hls_stream.h"

#define PRAGMA_SUB(PRAG) _Pragma(#PRAG)
#define CTX_PRAGMA(PRAG) PRAGMA_SUB(PRAG)

#ifndef CUR_SYN
    #define CUR_SYN 2
#endif

#ifndef CUR_SIM
    #define CUR_SIM 3
#endif

#define COND(NUM) \
    ((__SYNTHESIS__) && ((CUR_SYN) == (NUM))) || (!(__SYNTHESIS__) && ((CUR_SIM) == (NUM)))

#if COND(2)
    #define NUM_TROT         4
    #define NUM_SPIN         4096
    #define PACKET_SIZE      64
    #define LOG2_PACKET_SIZE 6
    #define NUM_STREAM       2
    #define LOG2_NUM_STREAM  1
    #define NUM_FADD         64
#elif COND(3)
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
/* Quantum Monte-Carlo U50 */
void QuantumMonteCarloU50(
    /* Spins */
    spin_pack_t trotters[NUM_TROT][NUM_SPIN / PACKET_SIZE],
    /* Jcoup */
    const fp_pack_t Jcoup_0[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
    const fp_pack_t Jcoup_1[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
    /* Array of h */
    const fp_t h[NUM_SPIN],
    /* Thermal Related */
    const fp_t Jperp,
    /* Thermal Related */
    const fp_t Beta,
    /* Log(Random Number = [0, 1]) * NUM_TROT */
    const fp_t logRand[NUM_TROT][NUM_SPIN]);
}

/* Naive Version */
void QuantumMonteCarloBasic(const int nTrot, const int nSpin, spin_t trotters[NUM_TROT][NUM_SPIN],
                            const fp_t Jcoup[NUM_SPIN][NUM_SPIN], const fp_t h[NUM_SPIN],
                            const fp_t Jperp, const fp_t Beta,
                            const fp_t logRandNumber[NUM_TROT][NUM_SPIN]);

#endif