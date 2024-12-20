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

    // モードを設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "MODE? (Simulation:0, Theory:1)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> MODE;

    sim.setSymbol();

    switch(sim.NUMBER_OF_BIT) {
        // QPSK
        case 2:
            // Simulation
            if(MODE == 0) {
                // ファイル作成
                filenameSimulation = "QPSKSimulation_Fading.csv";
                ofsSimulation.open(filenameSimulation);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    sim.setNoiseSD(EbN0dB);      // SNセット
                    ber = sim.getBerSimulation_Selective();

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

                    // ファイル出力
                    ofsSimulation << EbN0dB << "," << ber << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsSimulation.close();
            } 
            // Theory
            else if(MODE == 1) {
                // ファイル作成
                filenameTheory = "QPSKTheory_Fading.csv";
                ofsTheory.open(filenameTheory);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    berTheory = sim.getQPSKTheory_fading(EbN0dB);

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "Theory : " << EbN0dB << "," << berTheory << std::endl;

                    // ファイル出力
                    ofsTheory << EbN0dB << "," << berTheory << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsTheory.close();
            } else {
                break;
            }
        break;
        // 16QAM
        case 4:
            // Simulation
            if(MODE == 0) {
                // ファイル作成
                filenameSimulation = "16QAMSimulation_Fading.csv";
                ofsSimulation.open(filenameSimulation);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    sim.setNoiseSD(EbN0dB);      // SNセット
                    ber = sim.getBerSimulation_Selective();

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

                    // ファイル出力
                    ofsSimulation << EbN0dB << "," << ber << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsSimulation.close();
            }
            // Theory
            else if(MODE == 1) {
                // ファイル作成
                filenameTheory = "16QAMTheory_Fading.csv";
                ofsTheory.open(filenameTheory);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    berTheory = sim.get16QAMTheory_fading(EbN0dB);

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "Theory : " << EbN0dB << "," << berTheory << std::endl;

                    // ファイル出力
                    ofsTheory << EbN0dB << "," << berTheory << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsTheory.close();
            } else {
                break;
            }
        break;
        // 64QAM
        case 6:
            // Simulation
            if(MODE == 0) {
                // ファイル作成
                filenameSimulation = "64QAMSimulation_Fading.csv";
                ofsSimulation.open(filenameSimulation);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    sim.setNoiseSD(EbN0dB);      // SNセット
                    ber = sim.getBerSimulation_Selective();

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

                    // ファイル出力
                    ofsSimulation << EbN0dB << "," << ber << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsSimulation.close();
            } 
            // Theory
            else if(MODE == 1) {
                // ファイル作成
                filenameTheory = "64QAMTheory_Fading.csv";
                ofsTheory.open(filenameTheory);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    berTheory = sim.get64QAMTheory_fading(EbN0dB);

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "Theory : " << EbN0dB << "," << berTheory << std::endl;

                    // ファイル出力
                    ofsTheory << EbN0dB << "," << berTheory << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsTheory.close();
            } else {
                break;
            }
        break;
        // 256QAM
        case 8:
            // Simulation
            if(MODE == 0) {
                // ファイル作成
                filenameSimulation = "256QAMSimulation_Fading.csv";
                ofsSimulation.open(filenameSimulation);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    sim.setNoiseSD(EbN0dB);      // SNセット
                    ber = sim.getBerSimulation_Selective();

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

                    // ファイル出力
                    ofsSimulation << EbN0dB << "," << ber << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsSimulation.close();
            } 
            // Theory
            else if(MODE == 1) {
                // ファイル作成
                filenameTheory = "256QAMTheory_Fading.csv";
                ofsTheory.open(filenameTheory);

                for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                    berTheory = sim.get256QAMTheory_fading(EbN0dB);

                    // 標準出力
                    std::cout << "--------------------------------------------" << std::endl;
                    std::cout << "Theory : " << EbN0dB << "," << berTheory << std::endl;

                    // ファイル出力
                    ofsTheory << EbN0dB << "," << berTheory << std::endl;
                }
                std::cout << "--------------------------------------------" << std::endl;
                ofsTheory.close();
            } else {
                break;
            }
        break;
    }

    return 0;
}
