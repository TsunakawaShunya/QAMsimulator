/*
 * File:   Simulator.h
 * Author: Tsunakawa
 *
 * Created on 2024/07/02, 12:45
*/

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES
#include <C:/eigen-3.4.0/Eigen/Dense>
#include <C:/eigen-3.4.0/Eigen/Eigen>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <bitset>
#include <functional>
#include "random_collection.h"

class Simulator {
    public:
    int NUMBER_OF_BIT;      // ビット数(QAM方式)

    // コンストラクタ
    Simulator() {
        // 次元を設定
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Number of Bit? (QPSK:2, 16QAM:4)" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cin >> NUMBER_OF_BIT;

        // 試行回数を設定
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "TRIAL?" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cin >> numberOfTrial_;

        numberOfSymbols_ = pow(2, NUMBER_OF_BIT);

        // リサイズ
        num_.resize(numberOfSymbols_);
        grayNum_.resize(numberOfSymbols_);
        symbol_.resize(numberOfSymbols_);
        phi_.resize(numberOfSymbols_);
        obj_.resize(numberOfSymbols_);

        // データ（0 ~ 2^M - 1）をセット
        setNum();
        
        // 乱数の設定 
        unitIntUniformRand_.init(0.0, numberOfSymbols_ - 1, seed_);
        unitCNormalRand_.init(0.0, M_SQRT1_2, seed_);
        unitNormalRand_.init(0.0, 0.5, seed_);
    }

    // デストラクタ
    virtual ~Simulator() {
    }

    // シンボル設計（一般化）
    void setSymbol() {
        int M = numberOfSymbols_;               // M = 2^NUMBER_OF_BIT (シンボル数)
        int sqrtM = sqrt(M);                    // 実部/虚部のレベル数 (例: 16QAMならsqrtM=4)
        double P = 1.0 / (2.0 * (M - 1) / 3.0); // 平均送信電力の正規化

        for (int i = 0; i < numberOfSymbols_; i++) {
            int realBits = grayNum_[i] & ((1 << (NUMBER_OF_BIT / 2)) - 1);                          // ビット列の後半（実部）
            int imagBits = (grayNum_[i] >> (NUMBER_OF_BIT / 2)) & ((1 << (NUMBER_OF_BIT / 2)) - 1); // ビット列の前半（虚部）

            symbol_(i).real(-(2 * realBits - (sqrtM - 1)) * sqrt(P));
            symbol_(i).imag(-(2 * imagBits - (sqrtM - 1)) * sqrt(P));
        }
    }

    // 加法性雑音の標準偏差を設定
    void setNoiseSD(double EbN0dB) {
        noiseSD_ = sqrt(pow(10.0, -0.1 * EbN0dB) / (double)NUMBER_OF_BIT);        // 入力された  Eb/N0 [dB] から変換
    }

    // -------------- 結果 --------------
    // シミュレーション
    double getBerSimulation_AWGN() {
        int berCnt = 0;     // ビット誤りの総数
        for(int tri = 0; tri < numberOfTrial_; tri++) {
            berCnt += getBitErrorCount_AWGN();
        }

        return (double)berCnt / (double)numberOfTrial_ / (double)NUMBER_OF_BIT;
    }

    double getBerSimulation_Flat() {
        int berCnt = 0;     // ビット誤りの総数
        for(int tri = 0; tri < numberOfTrial_; tri++) {
            berCnt += getBitErrorCount_Flat();
        }

        return (double)berCnt / (double)numberOfTrial_ / (double)NUMBER_OF_BIT;
    }

    double getBerSimulation_Selective() {
        int berCnt = 0;     // ビット誤りの総数
        for(int tri = 0; tri < numberOfTrial_; tri++) {
            berCnt += getBitErrorCount_Selective();
        }

        return (double)berCnt / (double)numberOfTrial_ / (double)NUMBER_OF_BIT;
    }


    // 理論値
    // QPSK理論値
    double getQPSK_AWGNTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);
        return Q(sqrt(2 * gamma_b));
    }

    // 16QAM理論値
    double get16QAM_AWGNTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);
        return (3.0 * Q(sqrt(4.0 / 5.0 * gamma_b)) 
                + 2.0 * Q(3.0 * sqrt(4.0 / 5.0 * gamma_b)) 
                - Q(5.0 * sqrt(4.0 / 5.0 * gamma_b))) / 4.0;
    }

    // QPSK理論値
    double getQPSKTheory_fading(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);
        double ber = (1 - sqrt(gamma_b / (1.0 + gamma_b))) / 2;
        return ber;
    }

    // 16QAM理論値
    double get16QAMTheory_fading(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);

        return 3.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 / 2.0 / gamma_b))
                + 1.0 / 4.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / gamma_b / 9.0))
                - 1.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / gamma_b / 25.0));
    }

    protected:
    int numberOfSymbols_;                               // シンボル数
    double P_;                                           // 1 / 平均送信電力

    // データ
    std::vector<int> num_;                              // 数値データベクトル
    std::vector<int> grayNum_;                          // グレイ符号ベクトル
    Eigen::VectorXd phi_;                               // 位相ベクトル
    Eigen::VectorXcd symbol_;                           // シンボルベクトル

    // 通信
    int txData_;                                        // 送信データ
    int rxData_;                                        // 受信データ
    double noiseSD_;                                    // 雑音の標準偏差
    int numberOfTrial_;                                 // 試行回数
    std::complex<double> y_;                            // 受信信号ベクトル
    std::complex<double> x_;                            // 送信信号ベクトル
    std::complex<double> h_;                            // 伝送路応答行列
    std::complex<double> n_;                            // 雑音ベクトル

    // 乱数用
    unsigned long int seed_ = 10;                       // seed値
    normal_distribution<> unitNormalRand_;              // 正規乱数
    uniform_int_distribution<> unitIntUniformRand_;     // int型一様乱数
    cnormal_distribution<> unitCNormalRand_;            // 複素正規乱数
    cnormal_distribution<> randomNoise_;                // 複素正規乱数

    // 距離
    Eigen::VectorXd obj_;                          // 最尤推定における最小化の目的関数

    // グレイ符号化
    int setGrayCode(int num) {
        num = num ^ (num >> 1);
        return num;
    }

    // 数値データとグレイ符号のデータのセット
    void setNum() {
        for (int i = 0; i < numberOfSymbols_; i++) {
            num_[i] = i;
            grayNum_[i] = setGrayCode(num_[i]);
        }
    }

    // 伝送路
    // AWGN伝送路の初期化
    void initAWGN() {
        h_ = 1.0;
    }

    // フラットフェージング伝送路の初期化
    void initFlat() {
        h_.real(unitNormalRand_());
        h_.imag(h_.real());
    }

    // 選択性フェージング伝送路の初期化
    void initSelective() {
        h_ = unitCNormalRand_();
    }

    // 送信信号
    void set_x() {
        txData_ = unitIntUniformRand_();        // 送信データを設定
        x_ = symbol_(txData_);                   // 送信信号を設定
    }

    // 雑音
    void set_n() {
        n_ = unitCNormalRand_();
    }

    // 受信
    void set_y() {
        y_ = x_ + noiseSD_ * n_ / h_;
    }

    // 最尤推定で復調
    void set_rxDataByML() {
        for(int i = 0; i < numberOfSymbols_; i++) {
            obj_(i) = std::norm(h_) * std::norm((y_ - symbol_(i)));
        }
        
        Eigen::VectorXd::Index minColumn;       // ノルムが最小な index（つまり受信データ）
        obj_.minCoeff(&minColumn);
        rxData_ = minColumn;
    }

    // ハミング距離計算
    int hammingDistance(int num1, int num2) {
        int ham = 0;
        int xorResult;
        int bitMask = 1;

        xorResult = num1 ^ num2;

        for(int i = 0; i < NUMBER_OF_BIT; i++) {
            ham += (xorResult & bitMask) >> i;
            bitMask <<= 1;
        }
        
        return ham;
    }

    // ビット誤り数をカウント
    int getBitErrorCount_AWGN() {
        set_x();                    // 送信信号ベクトル生成
        /* 伝送路の初期化
        *  AWGN 伝送路：initAWGN()
        *  フラットフェージング伝送路：initFlat()
        *  選択性フェージング伝送路：initSelective()
        */
        initAWGN();                 // 伝送路応答行列を用意
        set_n();                    // 雑音生成
        set_y();                    // 受信
        set_rxDataByML();           //復調

        // ハミング距離計算
        return hammingDistance(grayNum_[txData_], grayNum_[rxData_]);
    }

    int getBitErrorCount_Flat() {
        set_x();                    // 送信信号ベクトル生成
        /* 伝送路の初期化
        *  AWGN 伝送路：initAWGN()
        *  フラットフェージング伝送路：initFlat()
        *  選択性フェージング伝送路：initSelective()
        */
        initFlat();                 // 伝送路応答行列を用意
        set_n();                    // 雑音生成
        set_y();                    // 受信
        set_rxDataByML();           //復調

        // ハミング距離計算
        return hammingDistance(grayNum_[txData_], grayNum_[rxData_]);
    }

    int getBitErrorCount_Selective() {
        set_x();                    // 送信信号ベクトル生成
        /* 伝送路の初期化
        *  AWGN 伝送路：initAWGN()
        *  フラットフェージング伝送路：initFlat()
        *  選択性フェージング伝送路：initSelective()
        */
        initSelective();                 // 伝送路応答行列を用意
        set_n();                    // 雑音生成
        set_y();                    // 受信
        set_rxDataByML();           //復調

        // ハミング距離計算
        return hammingDistance(grayNum_[txData_], grayNum_[rxData_]);
    }


    // Q関数
    double Q(double x) {
        return 1.0 / 2.0 * std::erfc(x / sqrt(2.0));
    }
};
#endif /* SIMULATOR_H */