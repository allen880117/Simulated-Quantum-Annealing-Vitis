#include "sqa.hpp"

/* Quantum Monte-Carlo */
void RunSQASoftware(const int nTrot, const int nSpin,
                    qubit_t    trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                    const fp_t Jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                    const fp_t h[MAX_QUBIT_NUM], const fp_t jperp,
                    const fp_t Beta,
                    const fp_t logRandNumber[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
{
    /* Tunnel-related energy */
    // fp_t dHTunnel = 2.0f * jperp * nTrot;
    fp_t dHTunnel = 2.0f * 0.5f * jperp * (fp_t)nTrot;

/* Traverse all trotters */
LOOP_TROTTERS:
    for (int m = 0; m < MAX_TROTTER_NUM; m++) {
        /* Exit condition */
        if (m == nTrot) break;

    /* Traverse all spins in trotter m */
    LOOP_SPINS:
        for (int i = 0; i < MAX_QUBIT_NUM; i++) {
            /* Exit condition */
            if (i == nSpin) break;

            /* Get this spin */
            qubit_t this_spin = trotters[m][i];

            /* Local Field of spin(m, i, j) */
            fp_t dH = 0;

        /* Compute Energy from same trotter */
        LOOP_COUPLES:
            for (int j = 0; j < MAX_QUBIT_NUM; j++) {
                /* Exit Condition */
                if (j == nSpin) break;
                if (trotters[m][j])
                    dH += Jcoup[i][j];
                else
                    dH -= Jcoup[i][j];
            }

            /* Compute Engery from up and down trotter */
            int     up           = (m == 0) ? (nTrot - 1) : (m - 1);
            int     down         = (m == nTrot - 1) ? (0) : (m + 1);
            qubit_t up_trotter   = trotters[up][i];
            qubit_t down_trotter = trotters[down][i];
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

            /* Times -2 and itself */
            dH *= -2.0f;
            if (!this_spin) { dH = -dH; }

            /* Flip */
            if ((-Beta * dH) > logRandNumber[m][i]) {
                trotters[m][i] = (!trotters[m][i]);
            }
        }
    }
}
