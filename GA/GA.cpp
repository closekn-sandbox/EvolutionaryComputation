#include <bits/stdc++.h>
using namespace std;

#define M 20  // 個体数
#define D 5   // 解の次元

double Pm = 0.05;             // 突然変異確率5%
int Tmax = 100;               // 最大繰り返し回数

double Weight[D], Value[D];   // D次元の配列
double Wmax;                  // Max Weight

double X[M][D], Xnext[M][D];  // M個×D次元の配列
double F[M];                  // 評価関数値を格納
double Fbest = 0, Xbest[D];   // 軍全体の最適解
int p1, p2;                   // 親個体
double prob, r, sumf;         // ルーレット選択で用いる
int i, d, d1, d2, tmp;        // 2点交叉

int binary_rand();                  // バイナリ[0,1]をランダムに返却
int uniform_d(void);                // 0<=x<=D-1 整数乱数
double uniform(void);               // 0<=x<=1 一様乱数
double eval_func(double x[]);       // 評価関数
int roulette_selection(double f[]); // ルーレット選択処理

int main() {
  // 初期化
  int t;
  srand( (unsigned) time(NULL) );
  for ( int i = 0; i < M; i++ ) {
    for ( int d = 0; d < D; d++ ) {
      X[i][d] = binary_rand();  // バイナリ格納
    }
  }
  
  for ( t = 1; t <= Tmax; t++ ) {

    // 各個体の評価値F[i]を計算
    for ( int i = 0; i < M; i++ ) {
      F[i] = eval_func(X[i]);
      // 最良値Fbestと最適解Xbest[]の更新
      if ( F[i] > Fbest ) {
        Fbest = F[i];
        for ( int d = 0; d < D; d++ ) {
          Xbest[d] = X[i][d];
        }
      }
    }

    for ( int i = 0; i < M; i++ ) {

      // ルーレット選択で親個体p1,p2を選ぶ
      do  {
        p1 = roulette_selection(F);
        p2 = roulette_selection(F);
      } while ( p1 == p2 );

      // 交叉する次元d1,d2をランダムで生成
      do {
        d1 = uniform_d();
        d2 = uniform_d();
      } while ( d1 == d2 );

      if ( d1 > d2 ) {  // d1 < d2 となるように入れかえ
        tmp = d1; d1 = d2; d2 = tmp;
      }

      // 2点交叉により、子Xnext[i][]を生成
      for ( int d = 0; d < D; d++ ) {
        if ( d <= d1 || d > d2 ) {
          Xnext[i][d] = X[p1][d];
        } else {
          Xnext[i][d] = X[p2][d];
        }
      }
      // 突然変異処理
      for ( int d = 0; d < D; d++ ) {
        if ( rand() < Pm ) {  // 確率Pmで突然変異
          Xnext[i][d] = 1 - Xnext[i][d];
        }
      }
    }
    // XにXnextを上書きし、次週準備
    for ( int i = 0; i < M; i++ ) {
      for ( int d = 0; d < D; d++ ) {
        X[i][d] = Xnext[i][d];
      }
    }

  }

  // 出力
  cout << "解の目的関数値 Fg = " << Fbest << endl;
  cout << "最適解 Xbest = [";
  for ( int d = 0; d < D; d++ ) {
    cout << Xbest[d] << ' ';
  }
  cout << ']' << endl;
  
  return 0;
}

// 評価関数
double eval_func(double x[]) {
  double f = 0;
  for ( int d = 0; d < D; d++ ) {
    f += x[d];
  }
  return f;
}

// ルーレット選択処理
int roulette_selection(double f[]) {
  double csum[M]; // 累積和
  int pos;
  sumf = csum[0] = f[0];
  for ( int i = 1; i < M; i++ ) {
    sumf += f[i];
    csum[i] = f[i] + csum[i-1];
  }
  double r = (int) (rand() / (RAND_MAX / sumf));
  for ( pos = 0; pos < M-1; pos++ ) {
    if ( r <= csum[pos] ) { break; }
  }
  return pos;
}

// バイナリ[0,1]のランダムな返却
int binary_rand() {
  return (int) (rand() / (RAND_MAX / 2));
}

// 0<=x<=1 一様乱数
double uniform(void) {
  return (double)rand() / ((double)RAND_MAX+1.0);
}

// 0<=x<=D-1 整数乱数
int uniform_d(void) {
  return (int) (rand() / (RAND_MAX / D+1));
}