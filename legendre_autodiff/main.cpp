// main.cpp
#include "autodiff.hpp"
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

// Recursive Legendre polynomials using our simple AutoDiff
template <typename T>
void LegendrePolynomials(int n, const AutoDiff<T>& x, std::vector<AutoDiff<T>>& P)
{
    if (n < 0) {
        P.clear();
        return;
    }
    P.resize(n + 1);

    P[0] = AutoDiff<T>(T(1));           // P₀(x) = 1
    if (n == 0) return;

    P[1] = x;                           // P₁(x) = x

    for (int k = 2; k <= n; ++k) {
        // Bonnet’s recurrence: (2k-1)·x·P_{k-1} - (k-1)·P_{k-2}
        // P_k = -------------------------------------------------
        //                       k
        AutoDiff<T> term1 = AutoDiff<T>(T(2 * k - 1)) * x * P[k-1];
        AutoDiff<T> term2 = AutoDiff<T>(T(k - 1))     * P[k-2];
        P[k] = (term1 - term2) / AutoDiff<T>(T(k));
    }
}

int main()
{
    const int N_points   = 400;
    const int max_order  = 5;

    std::ofstream file("legendre.csv");
    if (!file.is_open()) {
        std::cerr << "Error: could not open legendre.csv for writing!\n";
        return 1;
    }

    file << std::fixed << std::setprecision(12);

    // Header
    file << "x";
    for (int k = 0; k <= max_order; ++k)
        file << ",P" << k << ",P" << k << "_prime";
    file << "\n";

    // Evaluate on [-1, 1]
    for (int i = 0; i <= N_points; ++i) {
        double x_val = -1.0 + 2.0 * i / N_points;  // from -1 to +1
        AutoDiff<double> x = AutoDiff<double>::variable(x_val);

        std::vector<AutoDiff<double>> P;
        LegendrePolynomials(max_order, x, P);

        file << x_val;
        for (int k = 0; k <= max_order; ++k) {
            file << "," << P[k].val
                 << "," << P[k].der;
        }
        file << "\n";
    }

    file.close();
    std::cout << "legendre.csv written successfully ("
              << N_points + 1 << " points)!\n";
    std::cout << "Now run: python3 plot_legendre.py\n";

    return 0;
}