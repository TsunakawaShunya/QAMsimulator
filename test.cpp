// SNに対するBER
#include "simulator.h"
#include <vector>

// SN比
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 30.0;       // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 1.0;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// ファイル名
std::string filenameSimulation;
std::string filenameTheory;

// ファイルの中身
std::ofstream ofsSimulation;
std::ofstream ofsTheory;

// BER
double ber;
double berTheory;

// シミュレーションか理論計算か（0:simulation, 1:theory）
int MODE;

int main() {
    Simulator sim;

    sim.setSymbol();
    sim.setSymbol_a();

    return 0;
}
