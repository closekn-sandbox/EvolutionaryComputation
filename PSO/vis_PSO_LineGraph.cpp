#include <bits/stdc++.h>
#include <unistd.h>
using namespace std;

#define GNUPLOT_PATH "/usr/local/bin/gnuplot"

#define M 30    // 粒子数
#define D 2     // 解の次元

double c = 1.494, w = 0.729;    // PSOのパラメタ
int Tmax = 1000;                // 最大繰り返し回数
double Cr = 1e-5;               // 終了条件
double Xmin = -5.0, Xmax = 5.0; // 範囲
double X[M][D];                 // 解配列
double V[M][D];                 // 速度配列
double F[M];                    // 評価関数値
double Fp[M], Xp[M][D];         // pbest
double Fg, Xg[D];               // gbest
double Fg1;

double uniform(void);           // 0<=x<=1 一様乱数
double uniform_x(void);         // Xmin<=x<=Xmax 一様乱数
double sphere(double x[]);      // sphere function
double rastrigin(double x[]);   // rastrigin function

int main() {
    // gnuplot
    FILE *gp, *fp;
    if ( (gp = popen(GNUPLOT_PATH, "w")) == NULL ) {
      cout << "gnuplot open error" << endl;
      exit(EXIT_FAILURE);
    }

    // 初期化
    srand( (unsigned) time(NULL) );
    Fg = DBL_MAX;
    for ( int d = 0; d < D; d++ ) { Xg[d] = uniform_x(); }
    for ( int m = 0; m < M; m++ ) {
        Fp[m] = DBL_MAX;
        for ( int d = 0; d < D; d++ ) {
            V[m][d] = 0;
            X[m][d] = uniform_x();
            Xp[m][d] = uniform_x();
        }
    }

    // PSO
    fp = fopen("move_point.dat", "w");
    int t;
    for ( t = 0; t < Tmax; t++ ) {
        for ( int i = 0; i < M; i++ ) {
          //F[i] = sphere(X[i]);
          F[i] = rastrigin(X[i]);
          if ( F[i] < Fp[i] ) {
            Fp[i] = F[i];
            for ( int d = 0; d < D; d++ ) {
              Xp[i][d] = X[i][d];
            }
            if ( Fp[i] < Fg ) {
              Fg = Fp[i];
              for ( int d = 0; d < D; d++ ) {
                Xg[d] = Xp[i][d];
              }
            }
          }
        }
        if ( t == 0 ) { Fg1 = Fg; }
        fprintf(fp, "%d %f\n", t, Fg);
        if ( Fg < Cr ) { break; }
        for ( int i = 0; i < M; i++ ) {
            for ( int d = 0; d < D; d++ ) {
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
    fclose(fp);

    // make graph
    fprintf(gp, "set xrange [%d:%d]\n", 0, t);
    fprintf(gp, "set yrange [%f:%f]\n", 0.0, Fg1+1.0);
    fprintf(gp, "unset key\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'LineGraph.png'\n");
    fprintf(gp, "set title 'y = Fg, x = t'\n");
    fprintf(gp, "plot 'move_point.dat' with linespoints\n");

    //
    fprintf(gp, "exit\n");
    fflush(gp);
    pclose(gp);

    // 出力
    cout << "終了時刻 t = " << t << endl;
    cout << "解の目的関数値 Fg = " << Fg << endl;
    cout << "解 Xg = [";
    for ( int d = 0; d < D; d++ ) { cout << " " << Xg[d]; }
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

// Xmin <= x <= Xmax
double uniform_x(void) {
    return 2*Xmax * ( (double)rand() / ((double)RAND_MAX+1.0) ) - Xmax;
}

// 0<=x<=1 一様乱数
double uniform(void) {
    return (double)rand() / ((double)RAND_MAX+1.0);
}