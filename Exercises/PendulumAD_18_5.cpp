#include <iostream>
#include <fstream> 
#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <cmath>
#include <autodiff.hpp>

using namespace ASC_ode;

class PendulumAD : public NonlinearFunction
{
private:
  double m_length;
  double m_gravity;

public:
  PendulumAD(double length, double gravity=9.81) : m_length(length), m_gravity(gravity) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    T_evaluate<double>(x, f);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    Vector<AutoDiff<2>> x_ad(2);
    Vector<AutoDiff<2>> f_ad(2);

    x_ad(0) = Variable<0>(x(0));
    x_ad(1) = Variable<1>(x(1));
    T_evaluate<AutoDiff<2>>(x_ad, f_ad);

    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
         df(i,j) = f_ad(i).deriv()[j];
  }

  template <typename T>
  void T_evaluate (VectorView<T> x, VectorView<T> f) const
  {
    f(0) = x(1);

    T factor = T(-m_gravity / m_length);
    f(1) = sin(x(0)) * factor;
  }
};


int main()
{
  double tend = 15.0;
  int steps = 1000;
  double tau = tend/steps;

  Vector<> y_exp = { M_PI+0.001, 0 };
  Vector<> y_crank = { M_PI+0.001, 0 };

  auto rhs = std::make_shared<PendulumAD>(1.0);



   ExplicitEuler explicit_stepper(rhs);
   CrankNicolson crank_stepper(rhs);

  std::ofstream explicit_outfile ("data/ExplicitExercise_18_5_mass_PendulumAD.txt");
  explicit_outfile << 0.0 << "  " << y_exp(0) << " " << y_exp(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     explicit_stepper.DoStep(tau, y_exp);

     explicit_outfile << (i+1) * tau << "  " << y_exp(0) << " " << y_exp(1) << std::endl;
  }

  std::ofstream crank_outfile ("data/CrankExercise_18_5_PendulumAD.txt");
  crank_outfile << 0.0 << "  " << y_crank(0) << " " << y_crank(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     crank_stepper.DoStep(tau, y_crank);

     crank_outfile << (i+1) * tau << "  " << y_crank(0) << " " << y_crank(1) << std::endl;
  }


}
