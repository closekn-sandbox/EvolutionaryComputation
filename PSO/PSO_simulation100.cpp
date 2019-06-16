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

double c = 1.494, w = 0.729;    // PSOのパラメタ
int Tmax = 1000;                // 最大繰り返し回数
double Cr = 1e-5;               // 終了条件
double Xmin = -5, Xmax = 5;     // 範囲
double X[M][D];                 // 解配列
double V[M][D];                 // 速度配列
double F[M];                    // 評価関数値
double Fp[M], Xp[M][D];         // pbest
double Fg, Xg[D];               // gbest

double PSO_sphere(int Dim);        // Sphere Function における PSO (引数 解の次元)
double PSO_rastrigin(int Dim);     // Rastrigin Function における PSO (引数 解の次元)

double uniform(void);           // 0<=x<=1 一様乱数
double uniform_x(void);         // Xmin<=x<=Xmax 一様乱数
double sphere(double x[]);      // Sphere Function
double rastrigin(double x[]);   // Sastrigin Function

double unbiased_dispersion(double tmp[], double ave);   // 不変分散 n = S
double dispersion(double tmp[], double ave);            // 分散 n = S
double average(double tmp[]);                           // tmp[]の平均 n = S
bool isdigit_str(string s);     // 文字列数値判定

//-- main
int main(int argc, char *argv[]) {
    // 引数不正
    if ( argc == 1 ) { cout << "argv[1] is nothing." << endl; exit(1); }
    string tmp = argv[1];
    if ( ! isdigit_str(tmp) ) { cout << "argv[1] is not int." << endl; exit(1); }
    // 宣言
    double Tave;    // 平均終了時間
    double Fg_s[S]; // 各シミュレーションにおけるFg
    double Fgave;   // 平均Fg
    int dim = stoi(tmp);

    //-- 2関数シミュレーション
    cout << "D = " << dim << endl;
    // Sphere func.
    Tave = 0.0;
    for ( int i = 0; i < S; i++ ) {
        Tave += PSO_sphere(dim);
        Fg_s[i] = Fg;
    }
    Tave /= S;
    Fgave = average(Fg_s);
    cout << "Sphere func." << endl;
    cout << " 平均Fg = " << Fgave << endl;
    cout << " 不変分散σ-2 = " << unbiased_dispersion(Fg_s, Fgave) << endl;
    cout << " 標本分散s2 = " << dispersion(Fg_s, Fgave) << endl;
    cout << " 平均終了時間Tend = " << Tave << endl;
    // Rastrigin func.
    Tave = 0.0;
    for ( int i = 0; i < S; i++ ) {
        Tave += PSO_rastrigin(dim);
        Fg_s[i] = Fg;
    }
    Tave /= S;
    Fgave = average(Fg_s);
    cout << "Rastrigin func." << endl;
    cout << " 平均Fg = " << Fgave << endl;
    cout << " 不変分散σ-2 = " << unbiased_dispersion(Fg_s, Fgave) << endl;
    cout << " 標本分散s2 = " << dispersion(Fg_s, Fgave) << endl;
    cout << " 平均終了時間Tend = " << Tave << endl;
    

    return 0;
}

//-- Sphere Function における PSO
double PSO_sphere(int Dim) {
    // 初期化
    srand( (unsigned) time(NULL) );
    Fg = DBL_MAX;
    for ( int d = 0; d < Dim; d++ ) { Xg[d] = uniform_x(); }
    for ( int m = 0; m < M; m++ ) {
        Fp[m] = DBL_MAX;
        for ( int d = 0; d < Dim; d++ ) {
            V[m][d] = 0;
            X[m][d] = uniform_x();
            Xp[m][d] = uniform_x();
        }
    }

    // PSO
    int t;
    for ( t = 0; t < Tmax; t++ ) {
        for ( int i = 0; i < M; i++ ) {
            F[i] = sphere(X[i]);
            if ( F[i] < Fp[i] ) {
                Fp[i] = F[i];
                for ( int d = 0; d < Dim; d++ ) {
                    Xp[i][d] = X[i][d];
                }
                if ( Fp[i] < Fg ) {
                    Fg = Fp[i];
                    for ( int d = 0; d < Dim; d++ ) {
                        Xg[d] = Xp[i][d];
                    }
                }
            }
        }
        if ( Fg < Cr ) { break; }
        for ( int i = 0; i < M; i++ ) {
            for ( int d = 0; d < Dim; d++ ) {
                double r1, r2;
                while ( 1 ) {
                    r1 = uniform();
                    r2 = uniform();
                    if ( r1 != r2 ) { break; }
                }
                V[i][d] = w*V[i][d] + c*r1*(Xp[i][d]-X[i][d]) + c*r2*(Xg[d]-X[i][d]);
                X[i][d] = X[i][d] + V[i][d];
            }
        }
    }

    // 返却
    return t;
}

//-- Rastrigin Function における PSO
double PSO_rastrigin(int Dim) {
    // 初期化
    srand( (unsigned) time(NULL) );
    Fg = DBL_MAX;
    for ( int d = 0; d < Dim; d++ ) { Xg[d] = uniform_x(); }
    for ( int m = 0; m < M; m++ ) {
        Fp[m] = DBL_MAX;
        for ( int d = 0; d < Dim; d++ ) {
            V[m][d] = 0;
            X[m][d] = uniform_x();
            Xp[m][d] = uniform_x();
        }
    }

    // PSO
    int t;
    for ( t = 0; t < Tmax; t++ ) {
        for ( int i = 0; i < M; i++ ) {
            F[i] = rastrigin(X[i]);
            if ( F[i] < Fp[i] ) {
                Fp[i] = F[i];
                for ( int d = 0; d < Dim; d++ ) {
                    Xp[i][d] = X[i][d];
                }
                if ( Fp[i] < Fg ) {
                    Fg = Fp[i];
                    for ( int d = 0; d < Dim; d++ ) {
                        Xg[d] = Xp[i][d];
                    }
                }
            }
        }
        if ( Fg < Cr ) { break; }
        for ( int i = 0; i < M; i++ ) {
            for ( int d = 0; d < Dim; d++ ) {
                double r1, r2;
                while ( 1 ) {
                    r1 = uniform();
                    r2 = uniform();
                    if ( r1 != r2 ) { break; }
                }
                V[i][d] = w*V[i][d] + c*r1*(Xp[i][d]-X[i][d]) + c*r2*(Xg[d]-X[i][d]);
                X[i][d] = X[i][d] + V[i][d];
            }
        }
    }

    // 返却
    return t;
}

// Sphere Function
double sphere(double x[]) {
    double sum = 0;
    for ( int i = 0; i < D; i++ ) {
        sum += pow(x[i], 2);
    }
    return sum;
}

// Rastrigin Function
double rastrigin(double x[]) {
    double sum = 0;
    for ( int i = 0; i < D; i++ ) {
        sum += pow(x[i], 2) - 10*cos(2*M_PI*x[i]) + 10;
    }
    return sum;
}

// Xmin <= x <= Xmax
double uniform_x(void) {
    return 2*Xmax * ( (double)rand() / ((double)RAND_MAX+1.0) ) - Xmax;
}

// 0<=x<=1 一様乱数
double uniform(void) {
    return (double)rand() / ((double)RAND_MAX+1.0);
}

// 不変分散 n = S
double unbiased_dispersion(double tmp[], double ave) {
    double udis = 0.0;
    for ( int i = 0; i < S; i++ ) {
        udis += pow(tmp[i]-ave, 2); 
    }
    return udis / (S-1);
}

// 分散 n = S
double dispersion(double tmp[], double ave) {
    double dis = 0.0;
    for ( int i = 0; i < S; i++ ) {
        dis += pow(tmp[i]-ave, 2);
    }
    return dis / S;
}

// tmp[] の平均 n = S
double average(double tmp[]) {
    double sum = 0.0;
    for ( int i = 0; i < S; i++ ) {
        sum += tmp[i];
    }
    return sum / S;
}

// 文字列数値判定
bool isdigit_str(string s) {
    int len = s.length();
    for ( int i = 0; i < len; i++ ) {
        if ( !isdigit(s[i]) ) { return false; }
    }
    return true;
}