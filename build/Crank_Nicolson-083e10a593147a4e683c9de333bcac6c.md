# Crank–Nicolson Method 
## Exercise 17.4.1

In this exercise we had to:

- implement the Crank–Nicolson method
- compare the three methods to each other
- use the three methods to solve an electronic network

## Crank–Nicolson

The Crank–Nicolson method has been implemented in  
[timestepper.hpp](/src/timestepper.hpp).

The method is defined as:

$$
y_{i+1} = y_i + \frac{\tau}{2}\big(f(t_i, y_i) + f(t_{i+1}, y_{i+1})\big)
$$

To achieve this, we adapted the code for the implicit and explicit Euler methods and implemented the formula above.

```cpp

  class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau_half;
    std::shared_ptr<ConstantFunction> m_yold;
    Vector<> m_vecf;
    std::shared_ptr<ConstantFunction> m_fold;
  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau_half(std::make_shared<Parameter>(0.0)) 
      , m_vecf(rhs->dimF())
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      m_fold = std::make_shared<ConstantFunction>(rhs->dimF());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = m_yold + m_tau_half * (m_fold + m_rhs) - ynew;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau_half->set(0.5 * tau);

      this->m_rhs->evaluate(y, m_vecf);
      m_fold->set(m_vecf);

      NewtonSolver(m_equ, y);
    }
  };
```

## Comparison
![](Plots/mass_spring_steps_50.png)
![](Plots/mass_spring_steps_100.png)
![](Plots/mass_spring_steps_500.png)

This is a simple mass–spring system, therefore the energy should be constant and the movement periodic. In state space, this would appear as a perfect ellipse. However, as we can see, this is not the case for the implicit and explicit Euler methods.  
The explicit method increases the speed at each step relatively quickly, and similarly, the implicit method decreases the velocity with each step, almost as if the system were damped. Both of these effects decrease with an increase in step size.

By far the best method in this example is the Crank–Nicolson method, as it seems to keep the energy in the system constant.

## RC-Circuit

The RC circuit is modelled by the formula:
$$
U_0(t)=cos(100\pi t)
$$

$$
U_C(t) + R C \frac{dU_C}{dt}(t) = U_0(t)
$$

To bring this into an autonomous form, we treat $t$ as a state variable and set

$$
x = t.
$$

Then the right-hand side becomes

$$
U_0(x) = \cos(100\pi x),
$$

and the ODE can be rewritten as

$$
U_C'(t) = \frac{\cos(100\pi x) - U_C(t)}{RC},
$$

$$
x'(t) = 1.
$$

This RC-Circuit was implemented as a class in [timestepper.hpp](/demos/Exercise_17_4_1_RCCircuit.cpp)

```cpp
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
```
### Plots
In the plots, one can observe a similar trend as in the mass–spring system.The explicit Euler method has too much energy, while the implicit Euler method has too little.
However, unlike in the mass–spring case, both numerical solutions remain periodic, so the energy no longer grows or decays over time.

On the left we set the parameters to R=1, C=1 and on the right R=100, C=10^-6

<p align="center">
  <img src="Plots/RC50stepsR1C1.png" width="45%">
  <img src="Plots/RC50stepsR100C-6.png" width="45%">
</p>
<p align="center">
  <img src="Plots/RC100stepsR1C1.png" width="45%">
  <img src="Plots/RC100stepsR100C-6.png" width="45%">
</p>
<p align="center">
  <img src="Plots/RC500stepsR1C1.png" width="45%">
  <img src="Plots/RC500stepsR100C-6.png" width="45%">
</p>




