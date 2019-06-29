#include <bits/stdc++.h>
using namespace std;

#define M 30  // 個体数
#define D 5   // 解の次元

double Cr = 0.9, Fw = 0.5;    // DEのパラメタ
int Tmax = 1000;              // 最大繰り返し回数
double Fend = 1e-5;           // 終了条件
double Xmin = -5, Xmax = 5;   // 範囲

double X[M][D], Xnew[M][D];   // M個×D次元の配列
double V[D], U[D];            // D次元
double F[M], Ftmp;            // 
double Fbest, Xbest[D];       // 

double uniform(void);           // 0<=x<=1 一様乱数
double uniform_x(void);         // Xmin<=x<=Xmax 一様乱数
int uniform_m(void);            // 0<=x<=M-1 整数乱数
int uniform_d(void);            // 0<=x<=D-1 整数乱数

double sphere(double x[]);      // sphere function
double rastrigin(double x[]);   // rastrigin function


int main() {
  // 初期化
  int t = 0;
  srand( (unsigned) time(NULL) );
  for ( int m = 0; m < M; m++ ) {
    for ( int d = 0; d < D; d++ ) {
      X[m][d] = uniform_x();  // Xiを初期化
    }
  }
  Fbest = DBL_MAX; // Fbestを最大

  for ( t = 1; t <= Tmax; t++ ) {
    for ( int i = 0; i < M; i++ ) {
      // F[i]
      F[i] = sphere(X[i]);
      //F[i] = rastrigin(X[i]);
      // 3個体選出
      int a, b, c;
      a = b = c = uniform_m();
      while ( a == i ) { a = uniform_m(); }
      while ( b == i || b == a ) { b = uniform_m(); }
      while ( c == i || c == a || c == b ) { c = uniform_m(); }
      // 突然変異 新しいV生成
      for ( int j = 0; j < D; j++ ) {
        V[j] = X[a][j] + Fw*(X[b][j] - X[c][j]);
      }
      // ランダムな次元を選出
      int Jr = uniform_d();
      // 交叉
      for ( int j = 0; j < D; j++ ) {
        double ri = uniform();
        if ( ri < Cr || j == Jr ) {
          U[j] = V[j];
        } else {
          U[j] = X[i][j];
        }
      }
      // F(U)の計算とFbestの更新
      Ftmp = sphere(U);
      //Ftmp = rastrigin(U);
      if ( Ftmp < F[i] ) {
        F[i] = Ftmp;
        for ( int j = 0; j < D; j++ ) { Xnew[i][j] = U[j]; }
        if ( F[i] < Fbest ) {
          Fbest = F[i];
          for ( int j = 0; j < D; j++ ) { Xbest[j] = Xnew[i][j]; }
        }
      } else {
        for ( int j = 0; j < D; j++ ) { Xnew[i][j] = X[i][j]; }
      }
    }
    // Xの更新
    for ( int i = 0; i < M; i++ ) {
      for ( int j = 0; j < D; j++ ) {
        X[i][j] = Xnew[i][j];
      }
    }
    // Fbestが終了条件を満たす
    if ( Fbest < Fend ) { break; }
  }

  // 出力
  cout << "終了時刻 t = " << t << endl;
  cout << "解の目的関数値 Fbest = " << Fbest << endl;
  cout << "解 Xbest = [";
  for ( int d = 0; d < D; d++ ) { cout << " " << Xbest[d]; }
  cout << "]" << endl;
  return 0;
}

// sphere function
double sphere(double x[]) {
    double sum = 0;
    for ( int i = 0; i < D; i++ ) {
        sum += pow(x[i], 2);
    }
    return sum;
}

// rastrigin function
double rastrigin(double x[]) {
    double sum = 0;
    for ( int i = 0; i < D; i++ ) {
        sum += pow(x[i], 2) - 10*cos(2*M_PI*x[i]) + 10;
    }
    return sum;
}

// 0<=x<=D-1 整数乱数
int uniform_d(void) {
  return (int) (rand() / (RAND_MAX / D+1));
}

// 0<=x<=M-1 整数乱数
int uniform_m(void) {
  return (int) (rand() / (RAND_MAX / M+1));
}

// Xmin <= x <= Xmax
double uniform_x(void) {
  return 2*Xmax * ( (double)rand() / ((double)RAND_MAX+1.0) ) - Xmax;
}

// 0<=x<=1 一様乱数
double uniform(void) {
  return (double)rand() / ((double)RAND_MAX+1.0);
}