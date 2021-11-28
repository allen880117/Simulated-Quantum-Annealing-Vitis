#include "../include/sqa.hpp"

/* General Inpput State for Run Final */
struct state_t {
    u32_t  i_pack;          // pack index of current spin
    u32_t  i_spin;          // spin index of current spin
    spin_t up_spin;         // spin from up trotter
    spin_t down_spin;       // spin from down trotter
    fp_t   h_local;         // cache h
    fp_t   log_rand_local;  // cache log rand
};

/* Fix Info for Run Final */
struct info_t {
    u32_t m;             // Number of this trotter
    fp_t  beta;          // beta
    fp_t  de_qefct;      // + qefct energy
    fp_t  neg_de_qefct;  // - qefct energy
};

/*
 * Multiply
 * - Spin (boolean) times Jcoup
 */
inline fp_t Multiply(spin_t spin, fp_t jcoup) {
#pragma HLS INLINE
    return ((!spin) ? (-(jcoup)) : (jcoup));
}

/*
 * Trotter Unit
 * - UpdateOfTrotters      : Sum up spin[j] * Jcoup[i][j]
 * - UpdateOfTrottersFinal : Add other terms and do the flip
 */
fp_t UpdateOfTrotters(const spin_pack_t trotters_local[NUM_SPIN / PACKET_SIZE],
                      const fp_pack_t   jcoup_local[NUM_SPIN / PACKET_SIZE]) {
    // Buffer for de
    fp_t de_tmp = 0.0f;

    // Sum up
SUM_UP:
    for (u32_t ofst = 0, pack_ofst = 0; ofst < NUM_SPIN / PACKET_SIZE / NUM_STREAM;
         ofst++, pack_ofst += NUM_STREAM) {
        // Pramgas: Pipeline and Confine the usage of fadd
        // CTX_PRAGMA(HLS ALLOCATION operation instances = fadd limit = 64)
        CTX_PRAGMA(HLS PIPELINE)

        // Unpack and Multiply
    UNPACK_STREAM:
        for (u32_t strm_ofst = 0; strm_ofst < NUM_STREAM; strm_ofst++) {
        UNPACK_PACK:
            for (u32_t spin_ofst = 0; spin_ofst < PACKET_SIZE; spin_ofst++) {
                // Multiply
                de_tmp += Multiply(trotters_local[pack_ofst + strm_ofst][spin_ofst],
                                   jcoup_local[pack_ofst + strm_ofst].data[spin_ofst]);
            }
        }
    }

    // Return
    return de_tmp;
}

void UpdateOfTrottersFinal(const u32_t stage, const info_t info, const state_t state, const fp_t de,
                           spin_pack_t trotters_local[NUM_SPIN / PACKET_SIZE]) {
    bool inside = (stage >= info.m && stage < NUM_SPIN + info.m);
    if (inside) {
        // Cache
        fp_t   de_tmp    = de;
        spin_t this_spin = trotters_local[state.i_pack][state.i_spin];

        // Add de_qefct
        bool same_dir = (state.up_spin == state.down_spin);
        if (same_dir) { de_tmp += (state.up_spin) ? info.neg_de_qefct : info.de_qefct; }

        // Times 2.0f then Add h_local
        de_tmp *= 2.0f;
        de_tmp += state.h_local;

        /*
         * Formula: - (-2) * spin(i) * deTmp > lrn / beta
         * EqualTo:          spin(i) * deTmp > lrn / Beta / 2
         */
        // Times this_spin
        if (!this_spin) { de_tmp = -(de_tmp); }

        // Flip and Return
        if ((de_tmp) > state.log_rand_local / info.beta * 0.5f) {
            trotters_local[state.i_pack][state.i_spin] = (~this_spin);
        }
    }
}

extern "C" {
void QuantumMonteCarloU50(spin_pack_t     trotters[NUM_TROT][NUM_SPIN / PACKET_SIZE],
                          const fp_pack_t jcoup[NUM_SPIN][NUM_SPIN / PACKET_SIZE],
                          const fp_t h[NUM_SPIN], const fp_t jperp, const fp_t beta,
                          const fp_t log_rand[NUM_TROT][NUM_SPIN]) {
    // Interface
#pragma HLS INTERFACE mode = m_axi bundle = gmem0 port = trotters
#pragma HLS INTERFACE mode = m_axi bundle = gmem1 port = jcoup
#pragma HLS INTERFACE mode = m_axi bundle = gmem2 port = h
#pragma HLS INTERFACE mode = m_axi bundle = gmem3 port = log_rand

    // Pragma: Aggreate for better throughput
#pragma HLS AGGREGATE compact = auto variable = jcoup

    // Local trotters
    spin_pack_t trotters_local[NUM_TROT][NUM_SPIN / PACKET_SIZE];
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = trotters_local
#pragma HLS ARRAY_PARTITION dim = 2 type = cyclic factor = 2 variable = trotters_local

    // Local jcoup
    fp_pack_t jcoup_local[NUM_TROT][NUM_SPIN / PACKET_SIZE];
#pragma HLS AGGREGATE compact = auto variable = jcoup_local
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = jcoup_local
#pragma HLS ARRAY_PARTITION dim = 2 type = cyclic factor = 2 variable = jcoup_local

    // input state and de
    state_t state[NUM_TROT];
    fp_t    de[NUM_TROT];
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = state
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = de

    // qefct-Related Energy
    const fp_t de_qefct     = jperp * ((fp_t)NUM_TROT);
    const fp_t neg_de_qefct = -(de_qefct);

    // Fix info of trotter units
    info_t info[NUM_TROT];
INIT_INFO:
    for (u32_t m = 0; m < NUM_TROT; m++) {
#pragma HLS UNROLL
        info[m].m            = m;
        info[m].beta         = beta;
        info[m].de_qefct     = de_qefct;
        info[m].neg_de_qefct = neg_de_qefct;
    }

    // Read trotters to local memory
READ_TROTTERS:
    for (u32_t m = 0; m < NUM_TROT; m++) {
    READ_TROTTERS_1:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE; ofst++) {
#pragma HLS PIPELINE
            trotters_local[m][ofst] = trotters[m][ofst];
        }
    }

    // Loop of stage
LOOP_STAGE:
    for (u32_t stage = 0; stage < (NUM_SPIN + NUM_TROT - 1); stage++) {

        // Update offset, h_local, log_rand_local
    UPDATE_INPUT_STATE:
        for (u32_t m = 0; m < NUM_TROT; m++) {
#pragma HLS UNROLL
            u32_t ofst = ((stage + NUM_SPIN - m) & (NUM_SPIN - 1));

            // offset
            state[m].i_pack = (ofst >> (LOG2_PACKET_SIZE));
            state[m].i_spin = (ofst & (PACKET_SIZE - 1));

            // h, lrn
            state[m].h_local        = h[ofst];
            state[m].log_rand_local = log_rand[m][ofst];

            // up/down spin
            u32_t up           = (m == 0) ? (NUM_TROT - 1) : (m - 1);
            u32_t down         = (m == NUM_TROT - 1) ? (0) : (m + 1);
            state[m].up_spin   = trotters_local[up][state[m].i_pack][state[m].i_spin];
            state[m].down_spin = trotters_local[down][state[m].i_pack][state[m].i_spin];
        }

        // Read New Jcuop[0]
    READ_JCOUP:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE; ofst++) {
#pragma HLS PIPELINE
            jcoup_local[0][ofst] = jcoup[stage & (NUM_SPIN - 1)][ofst];
        }

        // Run Trotter Units
    UPDATE_OF_TROTTERS:
        for (u32_t m = 0; m < NUM_TROT; m++) {
#pragma HLS UNROLL
            de[m] = UpdateOfTrotters(trotters_local[m], jcoup_local[m]);
        }

        // Run final step of Trotter Units
    UPDATE_OF_TROTTERS_FINAL:
        for (u32_t m = 0; m < NUM_TROT; m++) {
#pragma HLS UNROLL
            UpdateOfTrottersFinal(stage, info[m], state[m], de[m], trotters_local[m]);
        }

        // Shift down jcoup_local
    SHIFT_JCOUP:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE; ofst++) {
#pragma HLS PIPELINE
            for (i32_t m = NUM_TROT - 2; m >= 0; m--) {
#pragma HLS UNROLL
                jcoup_local[m + 1][ofst] = jcoup_local[m][ofst];
            }
        }
    }

    // Write trotters_local to host memory
WRITE_TROTTERS:
    for (u32_t m = 0; m < NUM_TROT; m++) {
    WRITE_TROTTERS_1:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE; ofst++) {
#pragma HLS PIPELINE
            trotters[m][ofst] = trotters_local[m][ofst];
        }
    }
}
}
