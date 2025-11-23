// exercise_17.2.2_solution.cpp
// Standalone C++ program that implements the mass-spring example
// It implements Explicit Euler and the Improved Euler (mid-point) method,
// runs experiments, writes CSV output files, and prints error/energy summaries.

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

struct State {
    double y; // position
    double v; // velocity
};

// Right-hand side for mass-spring with m=1, k=1 by default
inline State rhs(const State& s, double m = 1.0, double k = 1.0) {
    State r;
    r.y = s.v;
    r.v = -(k / m) * s.y;
    return r;
}

// Explicit Euler: y_{n+1} = y_n + tau * f(y_n)
std::vector<State> explicit_euler(const State& y0, double tau, int steps) {
    std::vector<State> traj(steps + 1);
    traj[0] = y0;
    State y = y0;
    for (int i = 0; i < steps; i++) {
        State f = rhs(y);
        y.y += tau * f.y;
        y.v += tau * f.v;
        traj[i + 1] = y;
    }
    return traj;
}

// Improved Euler (midpoint-like as given in exercise):
// y_tilde = y_n + tau/2 * f(y_n)
// y_{n+1} = y_n + tau * f(y_tilde)
std::vector<State> improved_euler(const State& y0, double tau, int steps) {
    std::vector<State> traj(steps + 1);
    traj[0] = y0;
    State y = y0;
    for (int i = 0; i < steps; i++) {
        State f = rhs(y);
        State ytilde;
        ytilde.y = y.y + 0.5 * tau * f.y;
        ytilde.v = y.v + 0.5 * tau * f.v;
        State f2 = rhs(ytilde);
        y.y += tau * f2.y;
        y.v += tau * f2.v;
        traj[i + 1] = y;
    }
    return traj;
}

// analytical solution for m=1,k=1 with initial (1,0): y = cos(t), v = -sin(t)
inline State analytical(double t) {
    State s;
    s.y = std::cos(t);
    s.v = -std::sin(t);
    return s;
}

std::vector<double> energy(const std::vector<State>& traj, double m = 1.0, double k = 1.0) {
    std::vector<double> E(traj.size());
    for (size_t i = 0; i < traj.size(); ++i)
        E[i] = 0.5 * m * traj[i].v * traj[i].v + 0.5 * k * traj[i].y * traj[i].y;
    return E;
}

// compute max-norm and RMS error for position and velocity comparing to analytical
void compute_errors(const std::vector<State>& traj, double tau, double& max_err_pos, double& max_err_vel, double& rms_pos, double& rms_vel) {
    max_err_pos = max_err_vel = 0.0;
    double sum2p = 0.0, sum2v = 0.0;
    int N = (int)traj.size();
    for (int i = 0; i < N; i++) {
        double t = i * tau;
        State s = analytical(t);
        double dp = std::abs(traj[i].y - s.y);
        double dv = std::abs(traj[i].v - s.v);
        if (dp > max_err_pos) max_err_pos = dp;
        if (dv > max_err_vel) max_err_vel = dv;
        sum2p += dp * dp;
        sum2v += dv * dv;
    }
    rms_pos = std::sqrt(sum2p / N);
    rms_vel = std::sqrt(sum2v / N);
}

int main(int argc, char** argv) {
    double T = 20.0; // default end time
    int steps = 200;  // default steps
    if (argc >= 2) T = std::stod(argv[1]);
    if (argc >= 3) steps = std::stoi(argv[2]);
    double tau = T / steps;

    // initial condition y(0)=1, v(0)=0
    State y0; y0.y = 1.0; y0.v = 0.0;

    std::cout << "Running mass-spring experiment: T=" << T << ", steps=" << steps << ", tau=" << tau << "\n";

    auto traj_exp = explicit_euler(y0, tau, steps);
    auto traj_imp = improved_euler(y0, tau, steps);

    // write CSV files: columns t,y,v
    std::ofstream fe("explicit_euler.csv");
    fe << std::setprecision(15);
    for (int i = 0; i <= steps; i++) {
        double t = i * tau;
        fe << t << "," << traj_exp[i].y << "," << traj_exp[i].v << "\n";
    }
    fe.close();

    std::ofstream fi("improved_euler.csv");
    fi << std::setprecision(15);
    for (int i = 0; i <= steps; i++) {
        double t = i * tau;
        fi << t << "," << traj_imp[i].y << "," << traj_imp[i].v << "\n";
    }
    fi.close();

    // energies
    auto Eexp = energy(traj_exp);
    auto Eimp = energy(traj_imp);
    std::ofstream feE("energy_comparison.csv");
    feE << std::setprecision(15);
    feE << "t, E_explicit, E_improved, E_exact\n";
    for (int i = 0; i <= steps; i++) {
        double t = i * tau;
        double Eexact = 0.5 * (std::sin(t) * std::sin(t) + std::cos(t) * std::cos(t)); // always 0.5 for (1,0) initial
        feE << t << "," << Eexp[i] << "," << Eimp[i] << "," << Eexact << "\n";
    }
    feE.close();

    // compute and print errors
    double max_ep, max_ev, rms_p, rms_v;
    compute_errors(traj_exp, tau, max_ep, max_ev, rms_p, rms_v);
    std::cout << "Explicit Euler errors: max|y-y*|=" << max_ep << ", max|v-v*|=" << max_ev
        << ", rms_y=" << rms_p << ", rms_v=" << rms_v << "\n";

    compute_errors(traj_imp, tau, max_ep, max_ev, rms_p, rms_v);
    std::cout << "Improved Euler errors: max|y-y*|=" << max_ep << ", max|v-v*|=" << max_ev
        << ", rms_y=" << rms_p << ", rms_v=" << rms_v << "\n";

    std::cout << "CSV output: explicit_euler.csv, improved_euler.csv, energy_comparison.csv\n";
    std::cout << "Plot time evolution and phase plots from CSV (e.g. using Python/matplotlib).\n";

    return 0;
}