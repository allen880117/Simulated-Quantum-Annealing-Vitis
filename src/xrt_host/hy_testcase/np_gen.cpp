#include <fstream>
#include <iostream>
#include <cmath>
#define NUM_TROT         256
#define NUM_SPIN         64
#define PACKET_SIZE      32
#define LOG2_PACKET_SIZE 5
#define NUM_STREAM       2
#define LOG2_NUM_STREAM  1
#define NUM_FADD         64

int main(int argc, char *argv[]) {
    float rand_num[NUM_SPIN];
    float jcoup[NUM_SPIN][NUM_SPIN];
    float h[NUM_SPIN];

    std::ofstream sches("hy_sche.txt");
    std::ofstream rns("hy_rn.txt");
    std::ofstream jcoups("hy_jcoup.txt");
    std::ofstream hs("hy_h.txt");

    for (int i = 0; i < NUM_SPIN; i++) {
        rand_num[i] = (float)(i + 1) / (float)32 / 4.0f;
        rns << rand_num[i] << "\n";
    }

    for (int i = 0; i < NUM_SPIN; i++) {
        for (int j = 0; j < NUM_SPIN; j++) {
            jcoup[i][j] = rand_num[i] * rand_num[j];
            jcoups << jcoup[i][j] << " ";
        }
        jcoups << "\n";
    }

    for (int i = 0; i < NUM_SPIN; i++) {
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
        float Jperp = -0.5 * T * log(tanh(gamma / (float)NUM_TROT / T));
        sches << Jperp << " " << beta << "\n";
    }
}