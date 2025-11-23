
#include <iostream>
#include <fstream> 
#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <cmath>

using namespace ASC_ode;


class RCCircuit: public NonlinearFunction
{
private:
  double R;
  double C;

public:
  RCCircuit(double R_, double C_) : R(R_), C(C_) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    double Uc = x(0);   // capacitor voltage
    double t  = x(1);   // time

    f(0) = (std::cos(100.0*t*M_PI) - Uc)/(R*C);
    f(1) = 1.0;
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    //derivatives with respect to the first state variable (Uc)
    df(0,0) = -1.0 / (R * C);
    df(1,0) = 0.0;
    //derivative with respect to the second state variable (t)
    df(0,1) = -100.0 * M_PI * std::sin(100.0 * x(1)) / (R * C);
    df(1,1) = 0.0;

  }
};



int main()
{
  double tend = 0.05;
  int steps = 50;
  double tau = tend/steps;

  Vector<> y_imp = { 0, 0 }; // initializer list;
  Vector<> y_exp = { 0, 0 };
  Vector<> y_crank = { 0, 0 };

  auto rhs = std::make_shared<RCCircuit>(100.0, 10^(-6));



   ImplicitEuler implicit_stepper(rhs);
   ExplicitEuler explicit_stepper(rhs);
   CrankNicolson crank_stepper(rhs);

  std::ofstream implicit_outfile ("data/ImplicitExercise_17_4_1_mass_RCCircuit.txt");
  implicit_outfile << 0.0 << "  " << y_imp(0) << " " << y_imp(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     implicit_stepper.DoStep(tau, y_imp);

     implicit_outfile << (i+1) * tau << "  " << y_imp(0) << " " << y_imp(1) << std::endl;
  }


  std::ofstream explicit_outfile ("data/ExplicitExercise_17_4_1_RCCircuit.txt");
  explicit_outfile << 0.0 << "  " << y_exp(0) << " " << y_exp(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     explicit_stepper.DoStep(tau, y_exp);

     explicit_outfile << (i+1) * tau << "  " << y_exp(0) << " " << y_exp(1) << std::endl;
  }

  std::ofstream crank_outfile ("data/CrankNicolsonExercise_17_4_1_RCCircuit.txt");
  crank_outfile << 0.0 << "  " << y_crank(0) << " " << y_crank(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     crank_stepper.DoStep(tau, y_crank);

     crank_outfile << (i+1) * tau << "  " << y_crank(0) << " " << y_crank(1) << std::endl;
  }
}
