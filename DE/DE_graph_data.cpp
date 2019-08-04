/*
    Sphere Function & Rastrigin Function 
    100回シミュレーション
    argv[1] は 解の次元数 intで指定
*/

#include <bits/stdc++.h>
using namespace std;

#define S 100   // シミュレーション回数
#define M 30    // 粒子数
#define D 20    // 解の次元(最大)
#define Tmax 1000 // 最大繰り返し回数

double Cr = 0.9, Fw = 0.5;    // DEのパラメタ
double Fend = 1e-5;           // 終了条件
double Xmin = -5, Xmax = 5;   // 範囲

double X[M][D], Xnew[M][D];   // M個×D次元の配列
double V[D], U[D];            // D次元
double F[M], Ftmp;            // 
double Fbest, Xbest[D];       //
double F_s[S][Tmax];          // [シミュレーション回数][t] の 最良解
int s;                     // 現シミュレーション回数

double uniform(void);           // 0<=x<=1 一様乱数
double uniform_x(void);         // Xmin<=x<=Xmax 一様乱数
int uniform_m(void);            // 0<=x<=M-1 整数乱数
int uniform_d(void);            // 0<=x<=D-1 整数乱数

double sphere(double x[]);      // sphere function
double rastrigin(double x[]);   // rastrigin function

double DE_sphere(int Dim);        // Sphere Function における DE (引数 解の次元)
double DE_rastrigin(int Dim);     // Rastrigin Function における DE (引数 解の次元)

double average(int t);                           // tmp[]の平均 n = S
bool isdigit_str(string s);         // 文字列数値判定
void init_arr(); // 配列初期化

int main(int argc, char *argv[]) {
  // 引数不正
  if ( argc == 1 ) { cout << "argv[1] is nothing." << endl; exit(1); }
  string tmp = argv[1];
  if ( ! isdigit_str(tmp) ) { cout << "argv[1] is not int." << endl; exit(1); }
  // 宣言
  double Tave;    // 平均終了時間
  double F_save[Tmax]; // 各tにおけるFの平均
  int dim = stoi(tmp);

  //-- 2関数シミュレーション
    cout << "D = " << dim << endl;
    // Sphere func.
    cout << "Sphere func" << endl;
    for ( int s = 0; s < S; s++ ) {
      init_arr();
    }
    Tave = 0.0;
    for ( s = 0; s < S; s++ ) {
        Tave += DE_sphere(dim);
    }
    Tave /= S;
    for ( int i = 0; i < (int)Tave; i++ ) {
      F_save[i] = average(i);
    }
    for ( int i = 0; i < (int)Tave; i++ ) {
      cout << i+1 << ' ' << F_save[i] << endl;
    }
    
    // Rastrigin func.
    cout << "Rastrigin func" << endl;
    for ( int s = 0; s < S; s++ ) {
      init_arr();
    }
    Tave = 0.0;
    for ( int s = 0; s < S; s++ ) {
        Tave += DE_rastrigin(dim);
    }
    Tave /= S;
    for ( int i = 0; i < (int)Tave; i++ ) {
      F_save[i] = average(i);
    }
    for ( int i = 0; i < (int)Tave; i++ ) {
      cout << i+1 << ' ' << F_save[i] << endl;
    }

    return 0;
}

//-- sphere
double DE_sphere(int Dim) {
  // 初期化
  int t = 0;
  srand( (unsigned) time(NULL) );
  for ( int m = 0; m < M; m++ ) {
    for ( int d = 0; d < Dim; d++ ) {
      X[m][d] = uniform_x();  // Xiを初期化
    }
  }
  Fbest = DBL_MAX; // Fbestを最大

  for ( t = 0; t < Tmax; t++ ) {
    for ( int i = 0; i < M; i++ ) {
      // F[i]
      F[i] = sphere(X[i]);
      // 3個体選出
      int a, b, c;
      a = b = c = uniform_m();
      while ( a == i ) { a = uniform_m(); }
      while ( b == i || b == a ) { b = uniform_m(); }
      while ( c == i || c == a || c == b ) { c = uniform_m(); }
      // 突然変異 新しいV生成
      for ( int j = 0; j < Dim; j++ ) {
        V[j] = X[a][j] + Fw*(X[b][j] - X[c][j]);
      }
      // ランダムな次元を選出
      int Jr = uniform_d();
      // 交叉
      for ( int j = 0; j < Dim; j++ ) {
        double ri = uniform();
        if ( ri < Cr || j == Jr ) {
          U[j] = V[j];
        } else {
          U[j] = X[i][j];
        }
      }
      // F(U)の計算とFbestの更新
      Ftmp = sphere(U);
      if ( Ftmp < F[i] ) {
        F[i] = Ftmp;
        for ( int j = 0; j < Dim; j++ ) { Xnew[i][j] = U[j]; }
        if ( F[i] < Fbest ) {
          Fbest = F[i];
          for ( int j = 0; j < Dim; j++ ) { Xbest[j] = Xnew[i][j]; }
        }
      } else {
        for ( int j = 0; j < Dim; j++ ) { Xnew[i][j] = X[i][j]; }
      }
    }
    // Xの更新
    for ( int i = 0; i < M; i++ ) {
      for ( int j = 0; j < Dim; j++ ) {
        X[i][j] = Xnew[i][j];
      }
    }
    F_s[s][t] = Fbest;
    // Fbestが終了条件を満たす
    if ( Fbest < Fend ) { break; }
  }

  return t;
}

//-- rastrigin
double DE_rastrigin(int Dim) {
  // 初期化
  int t = 0;
  srand( (unsigned) time(NULL) );
  for ( int m = 0; m < M; m++ ) {
    for ( int d = 0; d < Dim; d++ ) {
      X[m][d] = uniform_x();  // Xiを初期化
    }
  }
  Fbest = DBL_MAX; // Fbestを最大

  for ( t = 0; t < Tmax; t++ ) {
    for ( int i = 0; i < M; i++ ) {
      // F[i]
      F[i] = rastrigin(X[i]);
      // 3個体選出
      int a, b, c;
      a = b = c = uniform_m();
      while ( a == i ) { a = uniform_m(); }
      while ( b == i || b == a ) { b = uniform_m(); }
      while ( c == i || c == a || c == b ) { c = uniform_m(); }
      // 突然変異 新しいV生成
      for ( int j = 0; j < Dim; j++ ) {
        V[j] = X[a][j] + Fw*(X[b][j] - X[c][j]);
      }
      // ランダムな次元を選出
      int Jr = uniform_d();
      // 交叉
      for ( int j = 0; j < Dim; j++ ) {
        double ri = uniform();
        if ( ri < Cr || j == Jr ) {
          U[j] = V[j];
        } else {
          U[j] = X[i][j];
        }
      }
      // F(U)の計算とFbestの更新
      Ftmp = rastrigin(U);
      if ( Ftmp < F[i] ) {
        F[i] = Ftmp;
        for ( int j = 0; j < Dim; j++ ) { Xnew[i][j] = U[j]; }
        if ( F[i] < Fbest ) {
          Fbest = F[i];
          for ( int j = 0; j < Dim; j++ ) { Xbest[j] = Xnew[i][j]; }
        }
      } else {
        for ( int j = 0; j < Dim; j++ ) { Xnew[i][j] = X[i][j]; }
      }
    }
    // Xの更新
    for ( int i = 0; i < M; i++ ) {
      for ( int j = 0; j < Dim; j++ ) {
        X[i][j] = Xnew[i][j];
      }
    }
    F_s[s][t] = Fbest;
    // Fbestが終了条件を満たす
    if ( Fbest < Fend ) { break; }
  }

  return t;
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

// tmp[] の平均 n = S
double average(int t) {
    double sum = 0.0;
    int ct = 0;
    for ( int i = 0; i < S; i++ ) {
      if ( F_s[i][t] != -1 ) {
        ct++;
        sum += F_s[i][t];
      }
    }
    return sum / ct;
}

// 文字列数値判定
bool isdigit_str(string s) {
    int len = s.length();
    for ( int i = 0; i < len; i++ ) {
        if ( !isdigit(s[i]) ) { return false; }
    }
    return true;
}

// 配列初期化
void init_arr() {
  for ( int t = 0; t < Tmax; t++ ) {
    F_s[s][t] = -1;
  }
}