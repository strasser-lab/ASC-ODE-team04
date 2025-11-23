# Exercise in 17.4.1
In this exercise we had to:

- implement the Crank Nicolson method
- compare the three methods to each other
- use the three methods to solve a Electronic Network

## Crank Nicolson

The Crank Nicolson Method has been implemented in [timestepper.hpp](/src/timestepper.hpp)
The method is defined as:

$$
y_{i+1} = y_i + \frac{\tau}{2} \big( f(t_i, y_i) +  f(t_{i+1}, y_{i+1})\big)
$$

To achive this we have adapted the code for the Implicit and Explixit Euler Method and implemented the formula above.

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

![](docs/Plots/mass_spring_steps_50.png)
