#include "hls_math.h"
#include "hlslib/xilinx/Operators.h"
#include "hlslib/xilinx/TreeReduce.h"
#include "sqa.hpp"

struct state_t {
    u32_t   packIdx, qubitIdx;
    qubit_t upQubit, dwQubit;
    fp_t    h, rndNum;
};

struct info_t {
    u32_t trotIdx;
    fp_t  qEffectPos, qEffectNeg;
};

static fp_t negate(fp_t input)
{
#pragma HLS INLINE
    union {
        fp_t    fpData;
        uint8_t uintData[sizeof(fp_t)];
    } converter;
    converter.fpData = input;
    converter.uintData[sizeof(fp_t) - 1] ^= 0x80;
    return converter.fpData;
}

static fp_t multiply(qubit_t spin, fp_t jcoup)
{
#pragma HLS INLINE
    return ((!spin) ? (negate(jcoup)) : (jcoup));
}

static fp_t calcEnergy(const qubitPack_t qubitsCache[NUM_COL_QUBIT_CACHE],
                       const fpPack_t    jcoupCacheBank0[NUM_COL_JCOUP_MEM_BANK],
                       const fpPack_t    jcoupCacheBank1[NUM_COL_JCOUP_MEM_BANK])
{
#pragma HLS INLINE OFF
    fp_t packLevelBuf[NUM_COL_JCOUP_MEM_BANK];
    fp_t qubitLevelBuf[JCOUP_BANK_NUM][PACKET_SIZE];

CALC_ENERGY:
    for (u32_t colIdx = 0, packIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK;
         colIdx++, packIdx += JCOUP_BANK_NUM) {
        CTX_PRAGMA(HLS ALLOCATION operation instances = fadd limit = NUM_FADD)
#pragma HLS PIPELINE
        for (u32_t qubitIdx = 0; qubitIdx < PACKET_SIZE; qubitIdx++) {
#pragma HLS UNROLL
            qubitLevelBuf[0][qubitIdx] =
                multiply(qubitsCache[packIdx + 0][qubitIdx], jcoupCacheBank0[colIdx][qubitIdx]);
            qubitLevelBuf[1][qubitIdx] =
                multiply(qubitsCache[packIdx + 1][qubitIdx], jcoupCacheBank1[colIdx][qubitIdx]);
        }
        packLevelBuf[colIdx] =
            hlslib::TreeReduce<fp_t, hlslib::op::Add<fp_t>, PACKET_SIZE>(qubitLevelBuf[0]) +
            hlslib::TreeReduce<fp_t, hlslib::op::Add<fp_t>, PACKET_SIZE>(qubitLevelBuf[1]);
    }
    return hlslib::TreeReduce<fp_t, hlslib::op::Add<fp_t>, NUM_COL_JCOUP_MEM_BANK>(packLevelBuf);
}

static void updateQubit(const u32_t stage, const info_t info, const state_t state,
                        const fp_t energy, qubitPack_t qubitsCache[NUM_COL_QUBIT_CACHE])
{
    if ((stage < info.trotIdx) || stage >= MAX_QUBIT_NUM + info.trotIdx) return;

    fp_t    eTmp     = energy;
    qubit_t thisSpin = qubitsCache[state.packIdx][state.qubitIdx];

    bool same_dir = (state.upQubit == state.dwQubit);
    if (same_dir) { eTmp += (state.upQubit) ? info.qEffectNeg : info.qEffectPos; }

    eTmp *= 2.0f;
    eTmp += state.h;

    /*
     * Formula: - (-2) * spin(i) * deTmp > lrn / beta
     * EqualTo:          spin(i) * deTmp > lrn / Beta / 2
     */
    if (!thisSpin) { eTmp = negate(eTmp); }

    // if ((de_tmp) > state.log_rand_local / info.beta * 0.5f) {}
    if (eTmp > state.rndNum) { qubitsCache[state.packIdx][state.qubitIdx] = (~thisSpin); }
}

static void prefetchJcoup(i32_t    stage,
                          fpPack_t jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                          fpPack_t jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                          fpPack_t jcoupPrefetch0[NUM_COL_JCOUP_MEM_BANK],
                          fpPack_t jcoupPrefetch1[NUM_COL_JCOUP_MEM_BANK])
{
#pragma HLS INLINE
PREFETCH_JCOUP:
    for (u32_t colIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK; colIdx++) {
#pragma HLS PIPELINE
        u32_t ofst             = (stage + 1) & (MAX_QUBIT_NUM - 1);
        jcoupPrefetch0[colIdx] = jcoupMemBank0[ofst][colIdx];
        jcoupPrefetch1[colIdx] = jcoupMemBank1[ofst][colIdx];
    }
}

static void prefetchRndNum(i32_t stage, fp_t rndNumCache[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE],
                           fp_t rndNumPrefetch[MAX_TROTTER_NUM])
{
#pragma HLS INLINE
PREFETCH_RNDNUM:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL factor = 4
        u32_t colIdx            = (((stage + 1) + MAX_QUBIT_NUM - trotIdx) & (MAX_QUBIT_NUM - 1));
        rndNumPrefetch[trotIdx] = rndNumCache[trotIdx][colIdx];
    }
}

static void updateState(u32_t stage, state_t state[MAX_TROTTER_NUM], fp_t hCache[MAX_QUBIT_NUM],
                        fp_t        rndNumPrefetch[MAX_TROTTER_NUM],
                        qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
UPDATE_INPUT_STATE:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL factor = 4
        u32_t colIdx    = ((stage + MAX_QUBIT_NUM - trotIdx) & (MAX_QUBIT_NUM - 1));
        u32_t upTrotIdx = (trotIdx == 0) ? (MAX_TROTTER_NUM - 1) : (trotIdx - 1);
        u32_t dwTrotIdx = (trotIdx == MAX_TROTTER_NUM - 1) ? (0) : (trotIdx + 1);
        u32_t packIdx   = colIdx / PACKET_SIZE;
        u32_t qubitIdx  = colIdx % PACKET_SIZE;
        state[trotIdx]  = state_t{.packIdx  = packIdx,
                                 .qubitIdx = qubitIdx,
                                 .upQubit  = qubitsCache[upTrotIdx][packIdx][qubitIdx],
                                 .dwQubit  = qubitsCache[dwTrotIdx][packIdx][qubitIdx],
                                 .h        = hCache[colIdx],
                                 .rndNum   = rndNumPrefetch[trotIdx]};
    }
}

static void shiftJcoupCache(u32_t    stage,
                            fpPack_t jcoupCacheBank0[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK],
                            fpPack_t jcoupCacheBank1[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK],
                            fpPack_t jcoupPrefetch0[NUM_COL_JCOUP_MEM_BANK],
                            fpPack_t jcoupPrefetch1[NUM_COL_JCOUP_MEM_BANK])
{
SHIFT_DOWN_JCOUP_CACHE:
    for (u32_t colIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK; colIdx++) {
#pragma HLS PIPELINE
        for (i32_t trotIdx = MAX_TROTTER_NUM - 2; trotIdx >= 0; trotIdx--) {
#pragma HLS UNROLL
            jcoupCacheBank0[trotIdx + 1][colIdx] = jcoupCacheBank0[trotIdx][colIdx];
            jcoupCacheBank1[trotIdx + 1][colIdx] = jcoupCacheBank1[trotIdx][colIdx];
        }
        jcoupCacheBank0[0][colIdx] = jcoupPrefetch0[colIdx];
        jcoupCacheBank1[0][colIdx] = jcoupPrefetch1[colIdx];
    }
}

extern "C" {
void RunSQAHardwareOneStep(u32_t nTrotters, u32_t nQubits, fp_t jperp, fp_t hCache[MAX_QUBIT_NUM],
                           fp_t        rndNumCache[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE],
                           fpPack_t    jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                           fpPack_t    jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                           qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
#pragma HLS INLINE off

    // Pragma: Aggreate for better throughput
#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank0
#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank1
#pragma HLS AGGREGATE compact = auto variable = qubitsCache

    // Pragma: Partition
#if RESOURCE_CONSTRAINT
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = hCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = rndNumCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = qubitsCache
#else
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = hCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = qubitsCache
#endif

    // Jcoup Cache and Prefetch Storage
    fpPack_t jcoupCacheBank0[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK];
    fpPack_t jcoupCacheBank1[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK];
    fpPack_t jcoupPrefetch0[NUM_COL_JCOUP_MEM_BANK];
    fpPack_t jcoupPrefetch1[NUM_COL_JCOUP_MEM_BANK];
    fp_t     rndNumPrefetch[MAX_TROTTER_NUM];
#if RESOURCE_CONSTRAINT
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = jcoupCacheBank0
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = jcoupCacheBank1
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = rndNumPrefetch
#else
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = jcoupCacheBank0
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = jcoupCacheBank1
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumPrefetch
#endif

    // Quantum Effect and Input State and Info of PE
    fp_t    qEffectPos = jperp * ((fp_t)nTrotters);
    fp_t    qEffectNeg = negate(qEffectPos);
    state_t state[MAX_TROTTER_NUM];
    info_t  info[MAX_TROTTER_NUM];
#pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = state
#pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = info

INIT_INFO:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL factor = 4
        info[trotIdx] =
            info_t{.trotIdx = trotIdx, .qEffectPos = qEffectPos, .qEffectNeg = qEffectNeg};
    }

LOOP_STAGE:
    for (i32_t stage = -1; stage < (MAX_QUBIT_NUM + MAX_TROTTER_NUM - 1); stage++) {
#pragma HLS PIPELINE off
        if (stage != -1) updateState(stage, state, hCache, rndNumPrefetch, qubitsCache);
        prefetchRndNum(stage, rndNumCache, rndNumPrefetch);
        if (stage != -1)
            shiftJcoupCache(stage, jcoupCacheBank0, jcoupCacheBank1, jcoupPrefetch0,
                            jcoupPrefetch1);
        prefetchJcoup(stage, jcoupMemBank0, jcoupMemBank1, jcoupPrefetch0, jcoupPrefetch1);

        if (stage != -1) {
        UPDATE_QUBITS:
            for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
                fp_t e = calcEnergy(qubitsCache[trotIdx], jcoupCacheBank0[trotIdx],
                                    jcoupCacheBank1[trotIdx]);
                updateQubit(stage, info[trotIdx], state[trotIdx], e, qubitsCache[trotIdx]);
            }
        }
    }
}
}

static fp_t genRndNum(i32_t &seed)
{
#pragma HLS INLINE
    const int i4_huge = 2147483647;
    int       k;
    float     r;
    k    = seed / 127773;
    seed = 16807 * (seed - k * 127773) - k * 2836;
    if (seed < 0) { seed = seed + i4_huge; }
    r = (float)(seed)*4.656612875E-10;
    return r;
}

static void fillRndNumCache(fp_t  rndNumCache[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE],
                            i32_t seed[MAX_TROTTER_NUM], fp_t beta)
{
#if USING_STD_RNG
    static std::mt19937                  rng(12347);
    std::uniform_real_distribution<fp_t> unif(0.0, 1.0);
#endif

#pragma HLS INLINE off

FILL_RNDNUM_CACHE:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
        for (u32_t colIdx = 0; colIdx < NUM_COL_RNDNUM_CACHE; colIdx++) {
#if USING_STD_RNG
            fp_t tmp                     = unif(rng);
            rndNumCache[trotIdx][colIdx] = log(tmp) / beta / 2.0f;
#else
            rndNumCache[trotIdx][colIdx] = log(genRndNum(seed[trotIdx])) / beta / 2.0f;
#endif
        }
    }
}

#ifndef __SYNTHESIS__
qubitPack_t qubitsMemLogHW[MAX_STEP_NUM][MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE];
#endif

extern "C" {
void RunSQAHardware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                    fpPack_t jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fpPack_t jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fp_t hMem[MAX_QUBIT_NUM], fp_t jperpMem[MAX_STEP_NUM],
                    qubitPack_t qubitsMem[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
#pragma HLS INTERFACE mode = m_axi bundle = gmem0 port = jcoupMemBank0
#pragma HLS INTERFACE mode = m_axi bundle = gmem1 port = jcoupMemBank1
#pragma HLS INTERFACE mode = m_axi bundle = gmem2 port = hMem
#pragma HLS INTERFACE mode = m_axi bundle = gmem3 port = jperpMem
#pragma HLS INTERFACE mode = m_axi bundle = gmem4 port = qubitsMem

#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank0
#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank1
#pragma HLS AGGREGATE compact = auto variable = qubitsMem

    /* Caches and Static Memory */
    qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE];
    fp_t        hCache[MAX_QUBIT_NUM];
    fp_t        rndNumCache0[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE];
    fp_t        rndNumCache1[MAX_TROTTER_NUM][NUM_COL_RNDNUM_CACHE];
    i32_t       seed[MAX_TROTTER_NUM];
#if RESOURCE_CONSTRAINT
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = qubitsCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = hCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = rndNumCache0
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = rndNumCache1
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = seed
    #pragma HLS AGGREGATE compact = auto variable = qubitsCache
#else
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = qubitsCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = cyclic factor = 4 variable = hCache
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumCache0
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumCache1
    #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = seed
    #pragma HLS AGGREGATE compact = auto variable = qubitsCache
#endif

CACHE_QUBITS:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
        for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
#pragma HLS PIPELINE
            qubitsCache[trotIdx][colIdx] = qubitsMem[trotIdx][colIdx];

CACHE_H:
    for (u32_t colIdx = 0; colIdx < MAX_QUBIT_NUM; colIdx++)
#pragma HLS PIPELINE
        hCache[colIdx] = hMem[colIdx];

INIT_SEED:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
#pragma HLS UNROLL
        seed[trotIdx] = initRndNumSeed + trotIdx;

    /* Prefill Random Number Memory */
    fillRndNumCache(rndNumCache0, seed, beta);

SQA_MAIN_LOOP:
    for (u32_t step = 0; step < MAX_STEP_NUM; step++) {
#pragma HLS PIPELINE off
        if (step < nSteps) {
            fp_t jperp = jperpMem[step];
            if (step & (0x01)) {
                RunSQAHardwareOneStep(nTrotters, nQubits, jperp, hCache, rndNumCache1,
                                      jcoupMemBank0, jcoupMemBank1, qubitsCache);
                fillRndNumCache(rndNumCache0, seed, beta);
            } else {
                RunSQAHardwareOneStep(nTrotters, nQubits, jperp, hCache, rndNumCache0,
                                      jcoupMemBank0, jcoupMemBank1, qubitsCache);
                fillRndNumCache(rndNumCache1, seed, beta);
            }
#ifndef __SYNTHESIS__
            if ((step + 1) % 20 == 0) std::cout << (step + 1) << " iterations done..." << std::endl;
            for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
                for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
                    qubitsMemLogHW[step][trotIdx][colIdx] = qubitsCache[trotIdx][colIdx];
#endif
        }
    }

WRITEBACK_CACHE:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
        for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
#pragma HLS PIPELINE
            qubitsMem[trotIdx][colIdx] = qubitsCache[trotIdx][colIdx];
}
}
