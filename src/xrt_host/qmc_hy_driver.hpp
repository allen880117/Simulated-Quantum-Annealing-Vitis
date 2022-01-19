#ifndef _QMC_HY_DRIVER_HPP_
#define _QMC_HY_DRIVER_HPP_

#include <random>
#include <string>

#include "experimental/xrt_bo.h"
#include "experimental/xrt_device.h"
#include "experimental/xrt_kernel.h"
#include "hls_math.h"

#define NUM_TROT         256
#define NUM_SPIN         64
#define PACKET_SIZE      32
#define LOG2_PACKET_SIZE 5
#define NUM_STREAM       2
#define LOG2_NUM_STREAM  1
#define NUM_FADD         64

typedef unsigned int u32_t;
typedef int          i32_t;

typedef float fp_t;
typedef struct {
    fp_t data[PACKET_SIZE];
} fp_pack_t;

typedef ap_uint<1>           spin_t;
typedef ap_uint<PACKET_SIZE> spin_pack_t;

class qmc_hy_driver {
   public:
    int n_trot;

    xrt::kernel krnl;
    xrt::run run;

    xrt::bo bo_trotters;
    xrt::bo bo_Jcoup_0;
    xrt::bo bo_Jcoup_1;
    xrt::bo bo_h;
    xrt::bo bo_logRand;

    spin_pack_t* bo_trotters_map;
    fp_pack_t*   bo_Jcoup_0_map;
    fp_pack_t*   bo_Jcoup_1_map;
    fp_t*        bo_h_map;
    fp_t*        bo_logRand_map;

    std::random_device rd;

   public:
    qmc_hy_driver(const char* xclbin_path, const int n_trot, float jcoup[NUM_SPIN][NUM_SPIN],
                  float h[NUM_SPIN]) {
        int device_index = 0;

        // Set Number of Trotters
        this->set_n_trot(n_trot);

        // Open device and Show Device Name
        xrt::device device = xrt::device(device_index);
        std::cout << "[INFO][O] -> Open Device with Index: " << device_index << std::endl;
        std::cout << "[INFO][-] -> Device Name: " << device.get_info<xrt::info::device::name>()
                  << std::endl;

        // Load the binary (xclbin)
        xrt::uuid uuid = device.load_xclbin(xclbin_path);
        std::cout << "[INFO][O] -> Load binary(.xclbin) : " << xclbin_path << std::endl;

        // Build the kernel
        this->krnl = xrt::kernel(device, uuid, "QuantumMonteCarloU50");

        // Allocate input buffer in global memory
        std::cout << "[INFO][-] -> Allocate Buffer in Global Memory" << std::endl;

        size_t trots_size_in_bytes = NUM_TROT * NUM_SPIN / PACKET_SIZE * sizeof(spin_pack_t);
        size_t Jcoup_size_in_bytes =
            NUM_SPIN * NUM_SPIN / PACKET_SIZE / NUM_STREAM * sizeof(fp_pack_t);
        size_t h_size_in_bytes       = NUM_SPIN * sizeof(fp_t);
        size_t logRand_size_in_bytes = NUM_TROT * NUM_SPIN * sizeof(fp_t);

        this->bo_trotters = xrt::bo(device, trots_size_in_bytes, krnl.group_id(1));
        this->bo_Jcoup_0  = xrt::bo(device, Jcoup_size_in_bytes, krnl.group_id(2));
        this->bo_Jcoup_1  = xrt::bo(device, Jcoup_size_in_bytes, krnl.group_id(3));
        this->bo_h        = xrt::bo(device, h_size_in_bytes, krnl.group_id(4));
        this->bo_logRand  = xrt::bo(device, logRand_size_in_bytes, krnl.group_id(7));

        // Map the contents of the buffer object into host memory
        std::cout << "[INFO][-] -> Map the buffer into the host memory" << std::endl;
        this->bo_trotters_map = bo_trotters.map<spin_pack_t*>();
        this->bo_Jcoup_0_map  = bo_Jcoup_0.map<fp_pack_t*>();
        this->bo_Jcoup_1_map  = bo_Jcoup_1.map<fp_pack_t*>();
        this->bo_h_map        = bo_h.map<fp_t*>();
        this->bo_logRand_map  = bo_logRand.map<fp_t*>();

        // Create the test data of trotters
        std::cout << "[INFO][-] -> Create initial trotters" << std::endl;
        // Do Nothing

        // Create Jcoup Buffer
        std::cout << "[INFO][-] -> Fill Jcoup buffer" << std::endl;
        this->fill_jcoup_buffer(jcoup);

        // Create h Buffer
        std::cout << "[INFO][-] -> Fill h buffer" << std::endl;
        this->fill_h_buffer(h);

        // Initial Sync
        this->initial_sync();
    };

    void run_qmc(float jperp, float beta) {
        this->fill_random_number();
        this->run = this->krnl(this->n_trot, this->bo_trotters, this->bo_Jcoup_0, this->bo_Jcoup_1,
                               this->bo_h, jperp, beta, this->bo_logRand);
        this->run.wait();
    }

    void set_n_trot(int n_trot) { this->n_trot = n_trot; }

    void get_trotters(bool trotters[NUM_TROT][NUM_SPIN]) {
        // Sync trotters back
        this->bo_trotters.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        for (int i = 0; i < NUM_TROT; i++) {
            for (int j = 0; j < NUM_SPIN / PACKET_SIZE; j++) {
                for (int k = 0; k < PACKET_SIZE; k++) {
                    trotters[i][j * PACKET_SIZE + k] =
                        this->bo_trotters_map[i * NUM_SPIN / PACKET_SIZE + j][k];
                }
            }
        }
    }

   private:
    void fill_jcoup_buffer(float jcoup[NUM_SPIN][NUM_SPIN]) {
        for (int i = 0; i < NUM_SPIN; i++) {
            for (int j = 0; j < NUM_SPIN / PACKET_SIZE / NUM_STREAM; j++) {
                fp_pack_t tmp_0;
                fp_pack_t tmp_1;

                for (int k = 0; k < PACKET_SIZE; k++) {
                    tmp_0.data[k] = jcoup[i][(j * NUM_STREAM + 0) * PACKET_SIZE + k];
                    tmp_1.data[k] = jcoup[i][(j * NUM_STREAM + 1) * PACKET_SIZE + k];
                }

                this->bo_Jcoup_0_map[i * NUM_SPIN / PACKET_SIZE / NUM_STREAM + j] = tmp_0;
                this->bo_Jcoup_1_map[i * NUM_SPIN / PACKET_SIZE / NUM_STREAM + j] = tmp_1;
            }
        }
    }

    void fill_h_buffer(float h[NUM_SPIN]) {
        for (int i = 0; i < NUM_SPIN; i++) { this->bo_h_map[i] = h[i]; }
    }

    void fill_random_number() {
        std::random_device                   rd;
        std::mt19937                         rng(rd());
        std::uniform_real_distribution<fp_t> unif(0.0, 1.0);

        // Read Log Random Number for Flipping
        for (int k = 0; k < NUM_TROT * NUM_SPIN; k++) {
            this->bo_logRand_map[k] = log(unif(rng)); // * (this->n_trot);
        }

        // Sync bo_logRand
        this->bo_logRand.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }

    void initial_sync() {
        // Sync Input Buffers (which won't change by host in the iterations)
        this->bo_trotters.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_Jcoup_0.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_Jcoup_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_h.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }
};

#endif