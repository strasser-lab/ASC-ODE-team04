
#include <iostream>
#include <fstream> 
#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -stiffness/mass*x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -stiffness/mass;
  }
};


int main()
{
  double tend = 4*M_PI;
  int steps = 500;
  double tau = tend/steps;

  Vector<> y_imp = { 1, 0 }; // initializer list;
  Vector<> y_exp = { 1, 0 };
  Vector<> y_crank = { 1, 0 };

  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);



   ImplicitEuler implicit_stepper(rhs);
   ExplicitEuler explicit_stepper(rhs);
   CrankNicolson crank_stepper(rhs);

  std::ofstream implicit_outfile ("data/ImplicitExercise_17_4_1_mass_spring.txt");
  implicit_outfile << 0.0 << "  " << y_imp(0) << " " << y_imp(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     implicit_stepper.DoStep(tau, y_imp);

     implicit_outfile << (i+1) * tau << "  " << y_imp(0) << " " << y_imp(1) << std::endl;
  }


  std::ofstream explicit_outfile ("data/ExplicitExercise_17_4_1_mass_spring.txt");
  explicit_outfile << 0.0 << "  " << y_exp(0) << " " << y_exp(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     explicit_stepper.DoStep(tau, y_exp);

     explicit_outfile << (i+1) * tau << "  " << y_exp(0) << " " << y_exp(1) << std::endl;
  }

  std::ofstream crank_outfile ("data/CrankNicolsonExercise_17_4_1_mass_spring.txt");
  crank_outfile << 0.0 << "  " << y_crank(0) << " " << y_crank(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
     crank_stepper.DoStep(tau, y_crank);

     crank_outfile << (i+1) * tau << "  " << y_crank(0) << " " << y_crank(1) << std::endl;
  }




}
