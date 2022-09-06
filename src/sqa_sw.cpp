#include "sqa.hpp"

/* Quantum Monte-Carlo */
void RunSQASoftware(u32_t nTrotters, u32_t nQubits, qubit_t qubits[MAX_TROTTER_NUM][MAX_QUBIT_NUM],
                    fp_t jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM], fp_t h[MAX_QUBIT_NUM], fp_t jperp,
                    fp_t beta, fp_t logRndNum[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
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

            /* Times -2 and itself */
            dH *= -2.0f;
            if (!this_spin) { dH = -dH; }

            /* Flip */
            if ((-beta * dH) > logRndNum[m][i]) { qubits[m][i] = (!qubits[m][i]); }
        }
    }
}
