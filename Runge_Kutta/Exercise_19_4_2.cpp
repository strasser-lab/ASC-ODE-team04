#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// typedef for clarity
using Vec = vector<double>;

// === ODE: y' = -y, y(0) = 1, exact solution y(t) = exp(-t) ===
double f(double t, double y) { return -y; }
double exact(double t) { return exp(-t); }

// === Explicit Runge-Kutta stepper ===
double RKStep(double y, double t, double h, const vector<vector<double>>& A,
              const vector<double>& b, const vector<double>& c) {
    int s = b.size();
    vector<double> k(s);
    for (int i = 0; i < s; ++i) {
        double ti = t + c[i]*h;
        double yi = y;
        for (int j = 0; j < i; ++j) yi += h * A[i][j] * k[j];
        k[i] = f(ti, yi);
    }
    double y_new = y;
    for (int i = 0; i < s; ++i) y_new += h * b[i] * k[i];
    return y_new;
}

// === Implicit RK for scalar y' = f(t,y) (simple fixed-point iteration) ===
double ImplicitRKStep(double y, double t, double h,
                      const vector<vector<double>>& A,
                      const vector<double>& b,
                      const vector<double>& c,
                      int max_iter=10, double tol=1e-12) {
    int s = b.size();
    vector<double> k(s, f(t,y)); // initial guess
    for (int iter=0; iter<max_iter; ++iter) {
        vector<double> k_new(s);
        double err = 0;
        for (int i=0; i<s; ++i) {
            double yi = y;
            for (int j=0; j<s; ++j) yi += h * A[i][j] * k[j];
            k_new[i] = f(t + c[i]*h, yi);
            err = max(err, fabs(k_new[i]-k[i]));
        }
        k = k_new;
        if (err < tol) break;
    }
    double y_new = y;
    for (int i=0; i<s; ++i) y_new += h * b[i] * k[i];
    return y_new;
}

int main() {
    double T = 1.0;
    int N = 10;
    double h = T/N;

    // === RK2 (Midpoint) ===
    vector<vector<double>> A_RK2 = {{0,0}, {0.5,0}};
    vector<double> b_RK2 = {0,1};
    vector<double> c_RK2 = {0,0.5};

    // === RK4 ===
    vector<vector<double>> A_RK4 = {{0,0,0,0},
                                    {0.5,0,0,0},
                                    {0,0.5,0,0},
                                    {0,0,1,0}};
    vector<double> b_RK4 = {1.0/6,1.0/3,1.0/3,1.0/6};
    vector<double> c_RK4 = {0,0.5,0.5,1};

    // === Gauss-Legendre 2-stage ===
    double sqrt3 = sqrt(3)/6;
    vector<vector<double>> A_GL2 = {{0.25, 0.25-sqrt3}, {0.25+sqrt3,0.25}};
    vector<double> b_GL2 = {0.5,0.5};
    vector<double> c_GL2 = {0.5-sqrt3,0.5+sqrt3};

    // === Radau IIA 2-stage ===
    vector<vector<double>> A_Radau2 = {{5.0/12,-1.0/12},{3.0/4,1.0/4}};
    vector<double> b_Radau2 = {3.0/4,1.0/4};
    vector<double> c_Radau2 = {1.0/3,1.0};

    // Initial condition
    double y_RK2 = 1.0, y_RK4 = 1.0, y_GL2 = 1.0, y_Radau2 = 1.0;

    cout << fixed << setprecision(8);
    cout << " t      RK2       error       RK4       error       GL2       error       Radau2    error\n";
    for (int i=0; i<=N; ++i) {
        double t = i*h;
        double exact_y = exact(t);

        double err_RK2 = fabs(y_RK2 - exact_y);
        double err_RK4 = fabs(y_RK4 - exact_y);
        double err_GL2 = fabs(y_GL2 - exact_y);
        double err_Radau2 = fabs(y_Radau2 - exact_y);

        cout << setw(4) << t
             << setw(10) << y_RK2 << setw(10) << err_RK2
             << setw(10) << y_RK4 << setw(10) << err_RK4
             << setw(10) << y_GL2 << setw(10) << err_GL2
             << setw(10) << y_Radau2 << setw(10) << err_Radau2
             << "\n";

        if (i<N) {
            y_RK2 = RKStep(y_RK2,t,h,A_RK2,b_RK2,c_RK2);
            y_RK4 = RKStep(y_RK4,t,h,A_RK4,b_RK4,c_RK4);
            y_GL2 = ImplicitRKStep(y_GL2,t,h,A_GL2,b_GL2,c_GL2);
            y_Radau2 = ImplicitRKStep(y_Radau2,t,h,A_Radau2,b_Radau2,c_Radau2);
        }
    }

    return 0;
}
