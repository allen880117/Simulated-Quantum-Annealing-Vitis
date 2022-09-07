#include "sqa.hpp"

/* Quantum Monte-Carlo */
void RunSQASoftwareOneStep(u32_t nTrotters, u32_t nQubits, fp_t jperp, fp_t h[MAX_QUBIT_NUM],
                           fp_t    rndNum[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                           fp_t    jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                           qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    /* Tunnel-related energy */
    // fp_t dHTunnel = 2.0f * jperp * nTrot;
    fp_t dHTunnel = 2.0f * 0.5f * jperp * (fp_t)nTrotters;

/* Traverse all trotters */
LOOP_TROTTERS:
    for (int m = 0; m < MAX_TROTTER_NUM; m++) {
        /* Exit condition */
        if (m == nTrotters) break;

    /* Traverse all spins in trotter m */
    LOOP_SPINS:
        for (int i = 0; i < MAX_QUBIT_NUM; i++) {
            /* Exit condition */
            if (i == nQubits) break;

            /* Get this spin */
            qubit_t this_spin = qubits[m][i];

            /* Local Field of spin(m, i, j) */
            fp_t dH = 0;

        /* Compute Energy from same trotter */
        LOOP_COUPLES:
            for (int j = 0; j < MAX_QUBIT_NUM; j++) {
                /* Exit Condition */
                if (j == nQubits) break;
                if (qubits[m][j])
                    dH += jcoup[i][j];
                else
                    dH -= jcoup[i][j];
            }

            /* Compute Engery from up and down trotter */
            int     up           = (m == 0) ? (nTrotters - 1) : (m - 1);
            int     down         = (m == nTrotters - 1) ? (0) : (m + 1);
            qubit_t up_trotter   = qubits[up][i];
            qubit_t down_trotter = qubits[down][i];
            if (up_trotter == down_trotter) {
                if (up_trotter)
                    dH -= dHTunnel;
                else
                    dH += dHTunnel;
            }

            // Times 2
            dH *= 2.0f;

            // Add h[i]
            dH += h[i];

#if 0
            /* Times -2 and itself */
            dH *= -2.0f;
            if (!this_spin) { dH = -dH; }

            /* Flip */
            if ((-beta * dH) > logRndNum[m][i]) { qubits[m][i] = (!qubits[m][i]); }
#else
            /* Times -2 and itself */
            if (!this_spin) { dH = -dH; }

            /* Flip */
            if (dH > rndNum[m][i]) { qubits[m][i] = (!qubits[m][i]); }
#endif
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
                          i32_t seed[MAX_TROTTER_NUM], fp_t beta)
{
#if USING_STD_RNG
    static std::mt19937                  rng(12347);
    std::uniform_real_distribution<fp_t> unif(0.0, 1.0);
#endif

#pragma HLS INLINE off
    for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++) {
#pragma HLS UNROLL
        for (u32_t colIdx = 0; colIdx < NUM_COL_RNDNUM_MEM; colIdx++) {
#if USING_STD_RNG
            fp_t tmp                   = unif(rng);
            rndNumMem[trotIdx][colIdx] = log(unif(rng)) / beta / 2.0f;
#else
            rndNumMem[trotIdx][colIdx] = log(genRndNum(seed[trotIdx])) / beta / 2.0f;
#endif
        }
    }
}

qubit_t qubitsLogSW[MAX_STEP_NUM][MAX_TROTTER_NUM][MAX_QUBIT_NUM];

void RunSQASoftware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                    fp_t jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t hMem[MAX_QUBIT_NUM],
                    fp_t jperpMem[MAX_STEP_NUM], qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    /* Caches and Static Memory */
    fp_t  rndNumMem0[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM];
    fp_t  rndNumMem1[MAX_TROTTER_NUM][NUM_COL_RNDNUM_MEM];
    i32_t seed[MAX_TROTTER_NUM];

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
                RunSQASoftwareOneStep(nTrotters, nQubits, jperp, hMem, rndNumMem1, jcoup, qubits);
                fillRndNumMem(rndNumMem0, seed, beta);
            } else {
                RunSQASoftwareOneStep(nTrotters, nQubits, jperp, hMem, rndNumMem0, jcoup, qubits);
                fillRndNumMem(rndNumMem1, seed, beta);
            }

            for (u32_t trotIdx = 0; trotIdx < MAX_TROTTER_NUM; trotIdx++)
                for (u32_t colIdx = 0; colIdx < MAX_QUBIT_NUM; colIdx++)
                    qubitsLogSW[step][trotIdx][colIdx] = qubits[trotIdx][colIdx];
        }
    }
}
