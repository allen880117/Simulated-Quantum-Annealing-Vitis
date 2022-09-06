#include <cmath>
#include <fstream>
#include <iostream>
#define MAX_TROTTER_NUM     256
#define MAX_QUBIT_NUM       64
#define PACKET_SIZE         32
#define LOG2_PACKET_SIZE    5
#define JCOUP_BANK_NUM      2
#define LOG2_JCOUP_BANK_NUM 1
#define NUM_FADD            64

int main(int argc, char *argv[])
{
    float rand_num[MAX_QUBIT_NUM];
    float jcoup[MAX_QUBIT_NUM][MAX_QUBIT_NUM];
    float h[MAX_QUBIT_NUM];

    std::ofstream sches("hy_sche.txt");
    std::ofstream rns("hy_rn.txt");
    std::ofstream jcoups("hy_jcoup.txt");
    std::ofstream hs("hy_h.txt");

    for (int i = 0; i < MAX_QUBIT_NUM; i++) {
        rand_num[i] = (float)(i + 1) / (float)32 / 4.0f;
        rns << rand_num[i] << "\n";
    }

    for (int i = 0; i < MAX_QUBIT_NUM; i++) {
        for (int j = 0; j < MAX_QUBIT_NUM; j++) {
            jcoup[i][j] = rand_num[i] * rand_num[j];
            jcoups << jcoup[i][j] << " ";
        }
        jcoups << "\n";
    }

    for (int i = 0; i < MAX_QUBIT_NUM; i++) {
        h[i] = 0.0f;
        hs << h[i] << "\n";
    }

    // Iter Arguments
    const int   iter        = 100;    // default 500
    const float gamma_start = 2.5f;   // default 3.0f
    const float T           = 0.05f;  // default 0.3f
    const float beta        = 1.0f / T;

    for (int i = 0; i < iter; i++) {
        float gamma = gamma_start * (1.0f - ((float)i / (float)iter));
        float jperp = -0.5 * T * log(tanh(gamma / (float)MAX_TROTTER_NUM / T));
        sches << jperp << " " << beta << "\n";
    }
}