/*
f1(x) = x1^2 + x2^2 -1 = 0
f2(x) = x1^3 - x2
初期値(x1, x2)=(1, 0.5) において
ニュートン法を用いてf(x)=0となるxを求めよ
途中解xkとf1(xk),f2(xk)の値も記せ
終了条件は f1(xk)<10^-6 AND f2(xk)<10^-6
*/

#include <bits/stdc++.h>
using namespace std;

#define         EPS (0.000001)
#define   f1(x1,x2) (x1*x1 + x2*x2 - 1)
#define   f2(x1,x2) (x1*x1*x1 - x2)
#define df11(x1,x2) (2*x1)
#define df12(x1,x2) (2*x2)
#define df21(x1,x2) (3*x1*x1)
#define df22(x1,x2) (-1)

void newton(double x1_0, double x2_0); // 2元ニュートン法

void main() {
    newton(1, 0.5);
    newton(-1, -0.5);
}

void newton(double x1_0, double x2_0) {
    double x1_new, x1_old, x2_new, x2_old;
    int ct = 1;
    x1_new = x1_0; x2_new = x2_0;

    cout << "x0=(" << x1_0 << ", " << x2_0 << ")" << endl;
    while ( f1(x1_new,x2_new) > EPS || f2(x1_new,x2_new) > EPS ) {
        x1_old = x1_new; x2_old = x2_new;
        x1_new = x1_old - (
            ( df22(x1_old,x2_old)*f1(x1_old,x2_old) + -1*df12(x1_old,x2_old)*f2(x1_old,x2_old) ) /
            ( df11(x1_old,x2_old)*df22(x1_old,x2_old) - df12(x1_old,x2_old)*df21(x1_old,x2_old) )
        );
        x2_new = x2_old - (
            ( -1*df21(x1_old,x2_old)*f1(x1_old,x2_old) + df11(x1_old,x2_old)*f2(x1_old,x2_old) ) /
            ( df11(x1_old,x2_old)*df22(x1_old,x2_old) - df12(x1_old,x2_old)*df21(x1_old,x2_old) )
        );
        cout << "x" << ct << "=(" << x1_new << ", " << x2_new << ") f(x" << ct << ")=(" << f1(x1_new,x2_new) << ", " << f2(x1_new,x2_new) << ")" << endl;
        ct++;
    }

    printf("x= (%f, %f)\n", x1_new, x2_new);
    cout << "x= (" << x1_new << ", " << x2_new << ")" << endl;
}
