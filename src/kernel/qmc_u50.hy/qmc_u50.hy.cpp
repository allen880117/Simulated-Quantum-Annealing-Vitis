#include "../include/sqa.hpp"
/* _QMC_U50_HY_CPP */
/* Special Edition for Hao-Yu Chung */

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
 * Negate
 * - Negate the sign of single-precision float-point
 * - To reduce the usage of LUT (replace the xor)
 */
inline float Negate(float input) {
#pragma HLS  INLINE
    union {
        float    fp_data;
        uint32_t int_data;
    } converter;

    converter.fp_data = input;

    ap_uint<32> tmp = converter.int_data;
    tmp[31]         = (~tmp[31]);

    converter.int_data = tmp;

    return converter.fp_data;
}

/*
 * Multiply
 * - Spin (boolean) times Jcoup
 */
inline fp_t Multiply(spin_t spin, fp_t jcoup) {
#pragma HLS INLINE
    return ((!spin) ? (Negate(jcoup)) : (jcoup));
}

/*
 * ReduceIntra (TOP)(BUF_SIZE = PACKET_SIZE)
 * - Recursion using template meta programming
 * - Reduce Intra-Buffer
 */
template <u32_t BUF_SIZE, u32_t GAP_SIZE>
inline void ReduceIntra(fp_t fp_buffer[BUF_SIZE]) {
#pragma HLS INLINE
    // Next call
    ReduceIntra<BUF_SIZE, GAP_SIZE / 2>(fp_buffer);

    // Reduce Intra
REDUCE_INTRA:
    for (u32_t i = 0; i < BUF_SIZE; i += GAP_SIZE) {
#pragma HLS UNROLL
        fp_buffer[i] += fp_buffer[i + GAP_SIZE / 2];
    }
}

/*
 * ReduceIntra (BOTTOM)(BUF_SIZE = PACKET_SIZE)
 */
template <>
inline void ReduceIntra<PACKET_SIZE, 1>(fp_t fp_buffer[PACKET_SIZE]) {
    ;
}

/*
 * ReduceIntra (BOTTOM)(BUF_SIZE = NUM_SPIN / PACKET_SIZE / NUM_STREAM)
 */
#if (NUM_SPIN / PACKET_SIZE / NUM_STREAM != PACKET_SIZE)
template <>
inline void ReduceIntra<NUM_SPIN / PACKET_SIZE / NUM_STREAM, 1>(
    fp_t fp_buffer[NUM_SPIN / PACKET_SIZE / NUM_STREAM]) {
    ;
}
#endif

/*
 * ReduceInter (TOP)
 * - Recursion using template meta programming
 * - Reduce Inter-Buffers
 */
template <u32_t N_STRM>
inline void ReduceInter(fp_t fp_buffer[NUM_STREAM][PACKET_SIZE]) {
#pragma HLS INLINE
    // Next call
    ReduceInter<N_STRM / 2>(fp_buffer);

    // Reduce Inter
REDUCE_INTER:
    for (u32_t i = 0; i < NUM_STREAM; i += N_STRM) {
#pragma HLS UNROLL
        fp_buffer[i][0] += fp_buffer[i + N_STRM / 2][0];
    }
}

/*
 * ReduceInter (BOTTOM)
 */
template <>
inline void ReduceInter<1>(fp_t fp_buffer[NUM_STREAM][PACKET_SIZE]) {
    ;
}

/*
 * Trotter Unit
 * - UpdateOfTrotters      : Sum up spin[j] * Jcoup[i][j]
 * - UpdateOfTrottersFinal : Add other terms and do the flip
 */
fp_t UpdateOfTrotters(
    const spin_pack_t trotters_local[NUM_SPIN / PACKET_SIZE],
    const fp_pack_t   jcoup_local[NUM_SPIN / PACKET_SIZE / NUM_STREAM][NUM_STREAM]) {
    // Buffer for de
    fp_t de_tmp[NUM_SPIN / PACKET_SIZE / NUM_STREAM];

    // Sum up
SUM_UP:
    for (u32_t ofst = 0, pack_ofst = 0; ofst < NUM_SPIN / PACKET_SIZE / NUM_STREAM;
         ofst++, pack_ofst += NUM_STREAM) {
        // Pramgas: Pipeline and Confine the usage of fadd
        CTX_PRAGMA(HLS ALLOCATION operation instances = fadd limit = NUM_FADD)
        CTX_PRAGMA(HLS PIPELINE)

        // Buffer for source of adder
        fp_t fp_buffer[NUM_STREAM][PACKET_SIZE];

        // Unpack and Multiply
    UNPACK_STREAM:
        for (u32_t strm_ofst = 0; strm_ofst < NUM_STREAM; strm_ofst++) {
#pragma HLS UNROLL

        UNPACK_PACK:
            for (u32_t spin_ofst = 0; spin_ofst < PACKET_SIZE; spin_ofst++) {
#pragma HLS UNROLL
                // Multiply
                fp_buffer[strm_ofst][spin_ofst] =
                    Multiply(trotters_local[pack_ofst + strm_ofst][spin_ofst],
                             jcoup_local[ofst][strm_ofst].data[spin_ofst]);
            }

            // Reduce inside each fp_buffer
            ReduceIntra<PACKET_SIZE, PACKET_SIZE>(fp_buffer[strm_ofst]);
        }

        // Reduce between different fp_buffer
        ReduceInter<NUM_STREAM>(fp_buffer);

        // Write into de_tmp buffer
        de_tmp[ofst] = fp_buffer[0][0];
    }

    // Reduce between de_tmp buffers
    ReduceIntra<NUM_SPIN / PACKET_SIZE / NUM_STREAM, NUM_SPIN / PACKET_SIZE / NUM_STREAM>(de_tmp);

    // Return
    return de_tmp[0];
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
        if (!this_spin) { de_tmp = Negate(de_tmp); }

        // Flip and Return
        if ((de_tmp) > state.log_rand_local / info.beta * 0.5f) {
            trotters_local[state.i_pack][state.i_spin] = (~this_spin);
        }
    }
}

#define ACT_TU 4

extern "C" {
void QuantumMonteCarloU50(const u32_t     n_trot,
                          spin_pack_t     trotters[NUM_TROT][NUM_SPIN / PACKET_SIZE],
                          const fp_pack_t jcoup_0[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
                          const fp_pack_t jcoup_1[NUM_SPIN][NUM_SPIN / PACKET_SIZE / NUM_STREAM],
                          const fp_t h[NUM_SPIN], const fp_t jperp, const fp_t beta,
                          const fp_t log_rand[NUM_TROT][NUM_SPIN]) {
    // Interface
#pragma HLS INTERFACE mode = m_axi bundle = gmem0 port = trotters
#pragma HLS INTERFACE mode = m_axi bundle = gmem1 port = jcoup_0
#pragma HLS INTERFACE mode = m_axi bundle = gmem2 port = jcoup_1
#pragma HLS INTERFACE mode = m_axi bundle = gmem3 port = h
#pragma HLS INTERFACE mode = m_axi bundle = gmem4 port = log_rand

    // Pragma: Aggreate for better throughput
#pragma HLS AGGREGATE compact = auto variable = jcoup_0
#pragma HLS AGGREGATE compact = auto variable = jcoup_1

    // Local trotters
    spin_pack_t trotters_local[NUM_TROT][NUM_SPIN / PACKET_SIZE];
    CTX_PRAGMA(HLS ARRAY_PARTITION dim = 1 type = cyclic factor = ACT_TU variable = trotters_local)
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = trotters_local

    // Local jcoup
    fp_pack_t jcoup_local[NUM_TROT][NUM_SPIN / PACKET_SIZE / NUM_STREAM][NUM_STREAM];
#pragma HLS AGGREGATE compact = auto variable = jcoup_local
    CTX_PRAGMA(HLS ARRAY_PARTITION dim = 1 type = cyclic factor = ACT_TU variable = jcoup_local)
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = jcoup_local

    // Prefetch jcoup, h, and log_rand
    fp_pack_t jcoup_prefetch[NUM_SPIN / PACKET_SIZE / NUM_STREAM][NUM_STREAM];
    fp_t      h_prefetch[NUM_TROT];
    fp_t      log_rand_prefetch[NUM_TROT];
    CTX_PRAGMA(HLS ARRAY_PARTITION dim = 1 type = cyclic factor = ACT_TU variable = h_prefetch)
    CTX_PRAGMA(HLS ARRAY_PARTITION dim = 1 type = cyclic factor = ACT_TU variable = log_rand_prefetch)
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = h_prefetch
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = log_rand_prefetch

    // input state and de and fix info of trotter units
    state_t state[NUM_TROT];
    fp_t    de[NUM_TROT];
    info_t  info[NUM_TROT];
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = state
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = de
// #pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = info

    // qefct-Related Energy
    const fp_t de_qefct     = jperp * ((fp_t)n_trot);
    const fp_t neg_de_qefct = Negate(de_qefct);

    // Initialize infos
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

    // Prefetch Jcoup before the loop of stages
PREFETCH_JCOUP:
    for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE / NUM_STREAM; ofst++) {
#pragma HLS PIPELINE
        jcoup_prefetch[ofst][0] = jcoup_0[0][ofst];
        jcoup_prefetch[ofst][1] = jcoup_1[0][ofst];
    }

    // Prefetch h and lr
    h_prefetch[0]        = h[0];
    log_rand_prefetch[0] = log_rand[0][0];

    // Loop of stage
LOOP_STAGE:
    for (u32_t stage = 0; stage < (NUM_SPIN + NUM_TROT - 1); stage++) {
#pragma HLS PIPELINE off
        // Update offset, h_local, log_rand_local
    UPDATE_INPUT_STATE:
        for (u32_t m = 0; m < NUM_TROT; m++) {
            CTX_PRAGMA(HLS UNROLL factor=ACT_TU)
            u32_t ofst = ((stage + NUM_SPIN - m) & (NUM_SPIN - 1));
            u32_t up   = (m == 0) ? (n_trot - 1) : (m - 1);
            u32_t down = (m == n_trot - 1) ? (0) : (m + 1);

            state[m].i_pack         = (ofst >> (LOG2_PACKET_SIZE));
            state[m].i_spin         = (ofst & (PACKET_SIZE - 1));
            state[m].up_spin        = trotters_local[up][state[m].i_pack][state[m].i_spin];
            state[m].down_spin      = trotters_local[down][state[m].i_pack][state[m].i_spin];
            state[m].h_local        = h_prefetch[m];
            state[m].log_rand_local = log_rand_prefetch[m];
        }

        // Read h and log_rand
    READ_H_LR:
        for (u32_t m = 0; m < NUM_TROT; m++) {
            CTX_PRAGMA(HLS UNROLL factor=ACT_TU)
            u32_t ofst           = (((stage + 1) + NUM_SPIN - m) & (NUM_SPIN - 1));
            h_prefetch[m]        = h[ofst];
            log_rand_prefetch[m] = log_rand[m][ofst];
        }

        // Shift down jcoup_local
    SHIFT_JCOUP:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE / NUM_STREAM; ofst++) {
            for (i32_t m = NUM_TROT - 2; m >= 0; m--) {
                CTX_PRAGMA(HLS UNROLL factor=ACT_TU)
                jcoup_local[m + 1][ofst][0] = jcoup_local[m][ofst][0];
                jcoup_local[m + 1][ofst][1] = jcoup_local[m][ofst][1];
            }
            jcoup_local[0][ofst][0] = jcoup_prefetch[ofst][0];
            jcoup_local[0][ofst][1] = jcoup_prefetch[ofst][1];
        }

        // Read New Jcuop[0]
    READ_JCOUP:
        for (u32_t ofst = 0; ofst < NUM_SPIN / PACKET_SIZE / NUM_STREAM; ofst++) {
#pragma HLS PIPELINE
            jcoup_prefetch[ofst][0] = jcoup_0[(stage + 1) & (NUM_SPIN - 1)][ofst];
            jcoup_prefetch[ofst][1] = jcoup_1[(stage + 1) & (NUM_SPIN - 1)][ofst];
        }

        // Run Trotter Units
    UPDATE_OF_TROTTERS:
        for (u32_t m = 0; m < NUM_TROT / ACT_TU; m++) {
            for (u32_t t = 0; t < ACT_TU ; t++){
#pragma HLS UNROLL
                de[m * ACT_TU + t] = UpdateOfTrotters(trotters_local[m * ACT_TU + t], jcoup_local[m * ACT_TU + t]);
                UpdateOfTrottersFinal(stage, info[m * ACT_TU + t], state[m * ACT_TU + t], 
                                      de[m * ACT_TU + t], trotters_local[m * ACT_TU + t]);
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
