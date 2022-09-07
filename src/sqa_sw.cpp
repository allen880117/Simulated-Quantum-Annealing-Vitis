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

void RunSQASoftware(u32_t nTrotters, u32_t nQubits, u32_t nSteps, fp_t beta, i32_t initRndNumSeed,
                    fp_t jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t hMem[MAX_QUBIT_NUM],
                    fp_t jperpMem[MAX_STEP_NUM], qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    ;
}