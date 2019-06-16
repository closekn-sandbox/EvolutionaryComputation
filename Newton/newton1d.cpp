/*
f(x) = x^3 +2x^2 -5x +6 がある
初期値x0=-3,1 の二通りの場合において
ニュートン法を用いてf(x)=0となるxを求めよ
*/

#include <bits/stdc++.h>
using namespace std;

#define   EPS (0.001)
#define  f(x) (x*x*x + 2*x*x - 5*x + 6)
#define df(x) (3*x*x + 4*x - 5)

void newton(double x0); // 1元ニュートン法

void main() {
    newton(-3);
    //newton(1); x0=1 X1=-1 x2=1 x3=-1 ...
}

void newton(double x0) {
    double x_new, x_old;
    int ct = 1;
    x_new = x0;

    cout << "x0= " << x0 << endl;
    while ( f(x_new) > EPS || -1*EPS > f(x_new) ) {
        x_old = x_new;
        x_new = x_old - ( f(x_old) / df(x_old) );
        cout << "x" << ct << "= " << x_new << " f(x" << ct << ")= " << f(x_new) << endl;
        ct++;
    }

    cout << "x= " << x_new << endl;
}