# AWGN伝送路シミュレータ

## 概要
こんなプログラムがあればすいすい進むのにと思って，優しいこの僕が作ってあげました．（感謝しろ:slightly_smiling_face: ）
このプロジェクトはAWGN（Additive White Gaussian Noise）伝送路でのBER（Bit Error Rate）のシミュレーション値と理論値を求めるためのコード．QPSKと16QAMだけやってみたので，他の256QAMや研究で発明した新しい設計では，それに合わせて書いて．実験の進め方として，**理論計算（気合で手計算）→シミュレーション**で理論計算が合っているか確認という手順になる．（理論計算は先生に縋っていました）

## ファイル構成
- `random_collection.h`: 藤井先生が作ってくれた正規乱数用のヘッダーファイル（ありがたく使ったけど自作もできると思う）
- `Simulator.h`: シミュレータクラスの宣言を含むヘッダーファイル
- `power_ber.cpp`: AWGN伝送路シミュレーションのメインプログラム

## 実行
プロジェクトのディレクトリに移ってこれをターミナルで実行
わからなければ，[vscode](https://code.visualstudio.com/)やらでC++の環境作って実行して
```
$ if ($?) { g++ power_ber.cpp -o power_ber } ; if ($?) { .\power_ber }
```

## ファイルの説明
### Simulator.h
* **呼び出されたらメンバ変数それぞれを設定**
* **シンボルをQPSK**
:::note warn
この時，平均のビット当たりのエネルギー(Eb)を1にしておく．
このコードではP_で制御してあります．
:::
* **雑音の標準偏差を設定**
これは大変気を付けた方がいいところです．1ビット当たりにかかる雑音のエネルギーをN0とすると，EbN0はEbは1になるようにしているから，N0はビットの数だけ割らないといけません．
* **シミュレーション**
* **理論値計算**

### power_ber.cpp
csvファイル作って結果を書く．