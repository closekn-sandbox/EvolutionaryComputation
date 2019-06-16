#include <bits/stdc++.h>

#define PI M_PI
#define LIM 100
#define f(x) (x*x*x + 2*x*x - 5*x + 6)

double uniform(void);                           // 0<=x<=1 一様乱数
double rand_normal(double mu, double sigma);    // 正規乱数(ガウス分布)
void markov_chain_algorithm();                  // マルコフチェインアルゴリズム

int main() {
    srand( (unsigned) time(NULL) );
    markov_chain_algorithm();
    return 0;
}

// マルコフチェインアルゴリズム
void markov_chain_algorithm() {
    double x_o, x_n;            // x_old, x_new
    double p,r;                 // %
    double optimal_solution;    // 最適解
    int t = 0;

    x_n = optimal_solution = uniform();
    printf("x%d = %f f(x%d) = %f\n", t, x_n, t, f(x_n));
    while ( fabs(f(x_n)) > 0.0001 && t <= LIM ) {
        t++;
        x_o = x_n;
        x_n = x_o + rand_normal(0,1);
        p = uniform();
        r = uniform();

        if ( r < p ) {
            // x_n = x_n;
        } else {
            x_n = x_o;
        }
        printf("x%d = %f f(x%d) = %f\n", t, x_n, t, f(x_n));
        if ( fabs(f(x_n)) < fabs(f(optimal_solution)) ) { optimal_solution = x_n; }
    }

    printf("Optimal Solution x = %f f(x) = %f\n", optimal_solution, f(optimal_solution));
}

// 0<=x<=1 一様乱数
double uniform(void) {
    return (double)rand() / ((double)RAND_MAX+1.0);
}

// 正規乱数(ガウス分布)
double rand_normal(double mu, double sigma) {
    double z = sqrt(-2.0*log(uniform())) * sin(2.0*PI*uniform());
    return mu + sigma*z;
}