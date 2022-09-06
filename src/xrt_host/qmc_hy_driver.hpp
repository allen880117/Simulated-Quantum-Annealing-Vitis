#ifndef _QMC_HY_DRIVER_HPP_
#define _QMC_HY_DRIVER_HPP_

#include <random>
#include <string>

#include "experimental/xrt_bo.h"
#include "experimental/xrt_device.h"
#include "experimental/xrt_kernel.h"
#include "hls_math.h"

#define MAX_TROTTER_NUM     256
#define MAX_QUBIT_NUM       64
#define PACKET_SIZE         32
#define LOG2_PACKET_SIZE    5
#define JCOUP_BANK_NUM      2
#define LOG2_JCOUP_BANK_NUM 1
#define NUM_FADD            64

typedef unsigned int u32_t;
typedef int          i32_t;

typedef float fp_t;
typedef struct {
    fp_t data[PACKET_SIZE];
} fpPack_t;

typedef ap_uint<1>           qubit_t;
typedef ap_uint<PACKET_SIZE> qubitPack_t;

class qmc_hy_driver
{
   public:
    int n_trot;

    xrt::kernel krnl;
    xrt::run    run;

    xrt::bo bo_trotters;
    xrt::bo bo_Jcoup_0;
    xrt::bo bo_Jcoup_1;
    xrt::bo bo_h;
    xrt::bo bo_logRand;

    qubitPack_t* bo_trotters_map;
    fpPack_t*    bo_Jcoup_0_map;
    fpPack_t*    bo_Jcoup_1_map;
    fp_t*        bo_h_map;
    fp_t*        bo_logRand_map;

    std::random_device rd;

   public:
    qmc_hy_driver(const char* xclbin_path, const int n_trot,
                  float jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM],
                  float h[MAX_QUBIT_NUM])
    {
        int device_index = 0;

        // Set Number of Trotters
        this->set_n_trot(n_trot);

        // Open device and Show Device Name
        xrt::device device = xrt::device(device_index);
        std::cout << "[INFO][O] -> Open Device with Index: " << device_index
                  << std::endl;
        std::cout << "[INFO][-] -> Device Name: "
                  << device.get_info<xrt::info::device::name>() << std::endl;

        // Load the binary (xclbin)
        xrt::uuid uuid = device.load_xclbin(xclbin_path);
        std::cout << "[INFO][O] -> Load binary(.xclbin) : " << xclbin_path
                  << std::endl;

        // Build the kernel
        this->krnl = xrt::kernel(device, uuid, "QuantumMonteCarloU50");

        // Allocate input buffer in global memory
        std::cout << "[INFO][-] -> Allocate Buffer in Global Memory"
                  << std::endl;

        size_t trots_size_in_bytes =
            MAX_TROTTER_NUM * MAX_QUBIT_NUM / PACKET_SIZE * sizeof(qubitPack_t);
        size_t Jcoup_size_in_bytes = MAX_QUBIT_NUM * MAX_QUBIT_NUM /
                                     PACKET_SIZE / JCOUP_BANK_NUM *
                                     sizeof(fpPack_t);
        size_t h_size_in_bytes = MAX_QUBIT_NUM * sizeof(fp_t);
        size_t logRand_size_in_bytes =
            MAX_TROTTER_NUM * MAX_QUBIT_NUM * sizeof(fp_t);

        this->bo_trotters =
            xrt::bo(device, trots_size_in_bytes, krnl.group_id(1));
        this->bo_Jcoup_0 =
            xrt::bo(device, Jcoup_size_in_bytes, krnl.group_id(2));
        this->bo_Jcoup_1 =
            xrt::bo(device, Jcoup_size_in_bytes, krnl.group_id(3));
        this->bo_h = xrt::bo(device, h_size_in_bytes, krnl.group_id(4));
        this->bo_logRand =
            xrt::bo(device, logRand_size_in_bytes, krnl.group_id(7));

        // Map the contents of the buffer object into host memory
        std::cout << "[INFO][-] -> Map the buffer into the host memory"
                  << std::endl;
        this->bo_trotters_map = bo_trotters.map<qubitPack_t*>();
        this->bo_Jcoup_0_map  = bo_Jcoup_0.map<fpPack_t*>();
        this->bo_Jcoup_1_map  = bo_Jcoup_1.map<fpPack_t*>();
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

    void run_qmc(float jperp, float beta)
    {
        this->fill_random_number();
        this->run = this->krnl(this->n_trot, this->bo_trotters,
                               this->bo_Jcoup_0, this->bo_Jcoup_1, this->bo_h,
                               jperp, beta, this->bo_logRand);
        this->run.wait();
    }

    void set_n_trot(int n_trot) { this->n_trot = n_trot; }

    void get_trotters(bool trotters[MAX_TROTTER_NUM][MAX_QUBIT_NUM])
    {
        // Sync trotters back
        this->bo_trotters.sync(XCL_BO_SYNC_BO_FROM_DEVICE);
        for (int i = 0; i < MAX_TROTTER_NUM; i++) {
            for (int j = 0; j < MAX_QUBIT_NUM / PACKET_SIZE; j++) {
                for (int k = 0; k < PACKET_SIZE; k++) {
                    trotters[i][j * PACKET_SIZE + k] =
                        this->bo_trotters_map[i * MAX_QUBIT_NUM / PACKET_SIZE +
                                              j][k];
                }
            }
        }
    }

   private:
    void fill_jcoup_buffer(float jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM])
    {
        for (int i = 0; i < MAX_QUBIT_NUM; i++) {
            for (int j = 0; j < MAX_QUBIT_NUM / PACKET_SIZE / JCOUP_BANK_NUM;
                 j++) {
                fpPack_t tmp_0;
                fpPack_t tmp_1;

                for (int k = 0; k < PACKET_SIZE; k++) {
                    tmp_0.data[k] =
                        jcoup[i][(j * JCOUP_BANK_NUM + 0) * PACKET_SIZE + k];
                    tmp_1.data[k] =
                        jcoup[i][(j * JCOUP_BANK_NUM + 1) * PACKET_SIZE + k];
                }

                this->bo_Jcoup_0_map[i * MAX_QUBIT_NUM / PACKET_SIZE /
                                         JCOUP_BANK_NUM +
                                     j] = tmp_0;
                this->bo_Jcoup_1_map[i * MAX_QUBIT_NUM / PACKET_SIZE /
                                         JCOUP_BANK_NUM +
                                     j] = tmp_1;
            }
        }
    }

    void fill_h_buffer(float h[MAX_QUBIT_NUM])
    {
        for (int i = 0; i < MAX_QUBIT_NUM; i++) { this->bo_h_map[i] = h[i]; }
    }

    void fill_random_number()
    {
        std::random_device                   rd;
        std::mt19937                         rng(rd());
        std::uniform_real_distribution<fp_t> unif(0.0, 1.0);

        // Read Log Random Number for Flipping
        for (int k = 0; k < MAX_TROTTER_NUM * MAX_QUBIT_NUM; k++) {
            this->bo_logRand_map[k] = log(unif(rng));  // * (this->n_trot);
        }

        // Sync bo_logRand
        this->bo_logRand.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }

    void initial_sync()
    {
        // Sync Input Buffers (which won't change by host in the iterations)
        this->bo_trotters.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_Jcoup_0.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_Jcoup_1.sync(XCL_BO_SYNC_BO_TO_DEVICE);
        this->bo_h.sync(XCL_BO_SYNC_BO_TO_DEVICE);
    }
};

#endif