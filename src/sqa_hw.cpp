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

template <u32_t STRIDE>
void reduceInterDim(fp_t buf[JCOUP_BANK_NUM][PACKET_SIZE])
{
    reduceInterDim<(STRIDE >> 1)>(buf);
    for (u32_t i = 0; i < JCOUP_BANK_NUM; i += STRIDE) {
#pragma HLS UNROLL
        buf[i][0] += buf[i + (STRIDE >> 1)][0];
    }
}

template <>
void reduceInterDim<1>(fp_t buf[JCOUP_BANK_NUM][PACKET_SIZE])
{
    ;
}

template <u32_t BUF_SIZE, u32_t STRIDE>
void reduceIntraDim(fp_t buf[BUF_SIZE])
{
#pragma HLS INLINE
    reduceIntraDim<BUF_SIZE, STRIDE / 2>(buf);
    for (u32_t i = 0; i < BUF_SIZE; i += STRIDE) {
#pragma HLS UNROLL
        buf[i] += buf[i + STRIDE / 2];
    }
}

template <>
void reduceIntraDim<PACKET_SIZE, 1>(fp_t buf[PACKET_SIZE])
{
    ;
}

#if (NUM_COL_JCOUP_MEM_BANK != PACKET_SIZE)
template <>
void reduceIntraDim<NUM_COL_JCOUP_MEM_BANK, 1>(fp_t buf[NUM_COL_JCOUP_MEM_BANK])
{
    ;
}
#endif

static fp_t calcEnergy(const qubitPack_t qubitsCache[NUM_COL_QUBIT_CACHE],
                       const fpPack_t    jcoupCache[NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM])
{
#pragma HLS INLINE OFF
    fp_t packLevelBuf[NUM_COL_JCOUP_MEM_BANK];
    fp_t qubitLevelBuf[JCOUP_BANK_NUM][PACKET_SIZE];

    for (u32_t colIdx = 0, packIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK;
         colIdx++, packIdx += JCOUP_BANK_NUM) {
        CTX_PRAGMA(HLS ALLOCATION operation instances = fadd limit = NUM_FADD)
#pragma HLS PIPELINE
        for (u32_t bankIdx = 0; bankIdx < JCOUP_BANK_NUM; bankIdx++) {
#pragma HLS UNROLL
            for (u32_t qubitIdx = 0; qubitIdx < PACKET_SIZE; qubitIdx++) {
#pragma HLS UNROLL
                qubitLevelBuf[bankIdx][qubitIdx] =
                    multiply(qubitsCache[packIdx + bankIdx][qubitIdx],
                             jcoupCache[colIdx][bankIdx].data[qubitIdx]);
            }
            reduceIntraDim<PACKET_SIZE, PACKET_SIZE>(qubitLevelBuf[bankIdx]);
        }
        reduceInterDim<JCOUP_BANK_NUM>(qubitLevelBuf);
        packLevelBuf[colIdx] = qubitLevelBuf[0][0];
    }
    reduceIntraDim<NUM_COL_JCOUP_MEM_BANK, NUM_COL_JCOUP_MEM_BANK>(packLevelBuf);
    return packLevelBuf[0];
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
                          fpPack_t jcoupPrefetch[NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM])
{
#pragma HLS INLINE
PREFETCH_JCOUP:
    for (u32_t colIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK; colIdx++) {
#pragma HLS PIPELINE
        u32_t ofst               = (stage + 1) & (MAX_QUBIT_NUM - 1);
        jcoupPrefetch[colIdx][0] = jcoupMemBank0[ofst][colIdx];
        jcoupPrefetch[colIdx][1] = jcoupMemBank1[ofst][colIdx];
    }
}

static void prefetchRndNum(i32_t stage, fp_t rndNumMem[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM],
                           fp_t rndNumPrefetch[MAX_TROTTER_NUM])
{
#pragma HLS INLINE
PREFETCH_RNDNUM:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
        u32_t colIdx            = (((stage + 1) + MAX_QUBIT_NUM - trotIdx) & (MAX_QUBIT_NUM - 1));
        rndNumPrefetch[trotIdx] = rndNumMem[trotIdx][colIdx];
    }
}

static void updateState(u32_t stage, state_t state[MAX_TROTTER_NUM], fp_t hCache[MAX_QUBIT_NUM],
                        fp_t        rndNumPrefetch[MAX_TROTTER_NUM],
                        qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
UPDATE_INPUT_STATE:
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
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

static void shiftJcoupCache(
    u32_t stage, fpPack_t jcoupCache[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM],
    fpPack_t jcoupPrefetch[NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM])
{
SHIFT_DOWN_JCOUP_CACHE:
    for (u32_t colIdx = 0; colIdx < NUM_COL_JCOUP_MEM_BANK; colIdx++) {
#pragma HLS PIPELINE
        for (i32_t trotIdx = MAX_TROTTER_NUM - 2; trotIdx >= 0; trotIdx--) {
#pragma HLS UNROLL
            jcoupCache[trotIdx + 1][colIdx][0] = jcoupCache[trotIdx][colIdx][0];
            jcoupCache[trotIdx + 1][colIdx][1] = jcoupCache[trotIdx][colIdx][1];
        }
        jcoupCache[0][colIdx][0] = jcoupPrefetch[colIdx][0];
        jcoupCache[0][colIdx][1] = jcoupPrefetch[colIdx][1];
    }
}

extern "C" {
void RunSQAHardware(u32_t nTrotters, u32_t nQubits, fp_t jperp, fp_t hCache[MAX_QUBIT_NUM],
                    fp_t        rndNumMem[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM],
                    fpPack_t    jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    fpPack_t    jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                    qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
#pragma HLS INLINE off
    // Pragma: Aggreate for better throughput
    // #pragma HLS AGGREGATE compact = auto variable = jcoupMemBank0
    // #pragma HLS AGGREGATE compact = auto variable = jcoupMemBank1
    // #pragma HLS AGGREGATE compact = auto variable = qubitsCache

    // Pragma: Partition
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumMem
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = qubitsCache

    // Jcoup Cache and Prefetch Storage
    fpPack_t jcoupCache[MAX_TROTTER_NUM][NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM];
    fpPack_t jcoupPrefetch[NUM_COL_JCOUP_MEM_BANK][JCOUP_BANK_NUM];
    fp_t     rndNumPrefetch[MAX_TROTTER_NUM];
// Disable these 2 pragmas, when this function is not top
// Don't know why these pragmas affect the memory partition and aggregation
// #pragma HLS AGGREGATE compact = auto variable = jcoupCache
// #pragma HLS AGGREGATE compact = auto variable = jcoupPrefetch
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = jcoupCache
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumPrefetch

    // Quantum Effect and Input State and Info of PE
    fp_t    qEffectPos = jperp * ((fp_t)nTrotters);
    fp_t    qEffectNeg = negate(qEffectPos);
    state_t state[MAX_TROTTER_NUM];
    info_t  info[MAX_TROTTER_NUM];
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = state
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = info

    // Initialize infos
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
        info[trotIdx] =
            info_t{.trotIdx = trotIdx, .qEffectPos = qEffectPos, .qEffectNeg = qEffectNeg};
    }

    // Prefetch Jcoup before the loop of stages
    prefetchJcoup(-1, jcoupMemBank0, jcoupMemBank1, jcoupPrefetch);
    prefetchRndNum(-1, rndNumMem, rndNumPrefetch);

LOOP_STAGE:
    for (u32_t stage = 0; stage < (MAX_QUBIT_NUM + MAX_TROTTER_NUM - 1); stage++) {
#pragma HLS PIPELINE off
        updateState(stage, state, hCache, rndNumPrefetch, qubitsCache);
        prefetchRndNum(stage, rndNumMem, rndNumPrefetch);
        shiftJcoupCache(stage, jcoupCache, jcoupPrefetch);
        prefetchJcoup(stage, jcoupMemBank0, jcoupMemBank1, jcoupPrefetch);

    UPDATE_QUBITS:
        for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
            fp_t e = calcEnergy(qubitsCache[trotIdx], jcoupCache[trotIdx]);
            updateQubit(stage, info[trotIdx], state[trotIdx], e, qubitsCache[trotIdx]);
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

static void fillRndNumMem(fp_t  rndNumMem[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM],
                          i32_t seed[MAX_TROTTER_NUM], u32_t beta)
{
#pragma HLS INLINE off
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
        for (u32_t colIdx = 0; colIdx < NUM_COL_RNDNUM_MEM; colIdx++) {
            rndNumMem[trotIdx][colIdx] = log(genRndNum(seed[trotIdx])) / beta / 2.0f;
        }
    }
}

#ifndef __SYNTHESIS__
qubitPack_t qubitsMemLogHW[MAX_STEP_NUM][MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE];
#endif

extern "C" {
void RunSQAKernel(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                  fpPack_t jcoupMemBank0[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                  fpPack_t jcoupMemBank1[MAX_QUBIT_NUM][NUM_COL_JCOUP_MEM_BANK],
                  fp_t hMem[MAX_QUBIT_NUM], fp_t jperpMem[MAX_STEP_NUM],
                  qubitPack_t qubitsMem[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE])
{
#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank0
#pragma HLS AGGREGATE compact = auto variable = jcoupMemBank1
#pragma HLS AGGREGATE compact = auto variable = qubitsMem

    /* Caches and Static Memory */
    qubitPack_t qubitsCache[MAX_TROTTER_NUM][NUM_COL_QUBIT_CACHE];
    fp_t        hCache[MAX_QUBIT_NUM];
    fp_t        rndNumMem0[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM];
    fp_t        rndNumMem1[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM];
    i32_t       seed[MAX_TROTTER_NUM];
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = qubitsCache
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumMem0
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = rndNumMem1
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = seed
#pragma HLS AGGREGATE compact = auto variable = qubitsCache

    /* Cache qubits */
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
        for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
#pragma HLS PIPELINE
            qubitsCache[trotIdx][colIdx] = qubitsMem[trotIdx][colIdx];

    /* Cache h */
    for (u32_t colIdx = 0; colIdx < MAX_QUBIT_NUM; colIdx++)
#pragma HLS PIPELINE
        hCache[colIdx] = hMem[colIdx];

    /* Init seed */
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
#pragma HLS UNROLL
        seed[trotIdx] = initRndNumSeed + trotIdx;

    /* Prefill Random Number Memory */
    fillRndNumMem(rndNumMem0, seed, beta);

    /* Main loop */
    for (u32_t step = 0; step < MAX_STEP_NUM; step++) {
#pragma HLS PIPELINE off
        if (step < nSteps) {
            fp_t jperp = jperpMem[step];
            if (step & (0x01)) {
                RunSQAHardware(nTrotters, nQubits, jperp, hCache, rndNumMem1, jcoupMemBank0,
                               jcoupMemBank1, qubitsCache);
                fillRndNumMem(rndNumMem0, seed, beta);
            } else {
                RunSQAHardware(nTrotters, nQubits, jperp, hCache, rndNumMem0, jcoupMemBank0,
                               jcoupMemBank1, qubitsCache);
                fillRndNumMem(rndNumMem1, seed, beta);
            }
#ifndef __SYNTHESIS__
            for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
                for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
                    qubitsMemLogHW[step][trotIdx][colIdx] = qubitsCache[trotIdx][colIdx];
#endif
        }
    }

    /* Writeback Cache */
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
        for (u32_t colIdx = 0; colIdx < NUM_COL_QUBIT_CACHE; colIdx++)
#pragma HLS PIPELINE
            qubitsMem[trotIdx][colIdx] = qubitsCache[trotIdx][colIdx];
}
}
