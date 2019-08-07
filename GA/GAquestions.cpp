#include <bits/stdc++.h>
using namespace std;

#define M 20  // 個体数
#define S 100 // シミュレーション回数
#define DIM 10  // 解の次元最大

int q_num;  // 問題番号
double Q1W[5] = { 7, 5, 1, 9, 6 };
double Q1V[5] = { 50, 40, 10, 70, 55 };
double Q2W[10] = { 3, 6, 5, 4, 8, 5, 3, 4, 8, 2 };
double Q2V[10] = { 70, 120, 90, 70, 130, 80, 40, 50, 30, 70 };
// シミュレーション内最適解
double allFb = 0;
double allXb[DIM];
int best_num = 0;

int D;          // 今回の解の次元数

double Pm = 0.0;             // 突然変異確率5%
int Tmax = 100;               // 最大繰り返し回数

double Weight[DIM], Value[DIM];   // D次元の配列
double Wmax;                  // Max Weight

double X[M][DIM], Xnext[M][DIM];  // M個×D次元の配列
double F[M];                  // 評価関数値を格納
double Fbest, Xbest[DIM];   // 軍全体の最適解
int p1, p2;                   // 親個体
double prob, r, sumf;         // ルーレット選択で用いる
int i, d, d1, d2, tmp;        // 2点交叉

bool setting_WV(int q);   // 今回の問題weight,valueのset
void GA();                // GA

int binary_rand();                  // バイナリ[0,1]をランダムに返却
int uniform_d(void);                // 0<=x<=D-1 整数乱数
double uniform(void);               // 0<=x<=1 一様乱数
double eval_func(double x[]);       // 評価関数
int roulette_selection(double f[]); // ルーレット選択処理

int main(int argc, char *argv[]) {
  // 初期化
  if ( argc < 2 ) { puts("arg[1] = the question number."); exit(1); }
  q_num = stoi(argv[1]);
  if ( ! setting_WV(q_num) ) { puts("There isn't the question."); exit(1); }
  srand( (unsigned) time(NULL) );

  // simulation
  for ( int i = 0; i < S; i++ ) {
    Fbest = 0;
    GA();
    if ( allFb == Fbest ) { best_num++; }
    if ( allFb < Fbest ) {
      allFb = Fbest;
      best_num = 1;
      for ( int d = 0; d < D; d++ ) {
        allXb[d] = Xbest[d];
      }
    }
  }

  // 出力
  cout << "Q" << q_num << endl;
  cout << "解の目的関数値 Fg = " << allFb << endl;
  cout << "最適解 Xbest = [";
  for ( int d = 0; d < D; d++ ) {
    cout << allXb[d] << ' ';
  }
  cout << ']' << endl;
  cout << "最適値を導き出した回数 : " << best_num << endl;
  
  return 0;
}

void GA() {
  int t;
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
}

bool setting_WV(int q) {
  switch ( q ) {
    case 1:
      D = 5;
      for ( int d = 0; d < D; d++ ) {
        Weight[d] = Q1W[d];
        Value[d] = Q1V[d];
      }
      Wmax = 15;
      return true;
    case 2:
      D = 10;
      for ( int d = 0; d < D; d++ ) {
        Weight[d] = Q2W[d];
        Value[d] = Q2V[d];
      }
      Wmax = 20;
      return true;
    default:
      return false;
  }
}

// 評価関数
double eval_func(double x[]) {
  double f = 0;
  double w = 0;
  for ( int d = 0; d < D; d++ ) {
    f += x[d]*Value[d];
    w += x[d]*Weight[d];
  }
  return ( w <= Wmax ) ? f : 1;
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