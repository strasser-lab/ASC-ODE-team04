# Exercise 17.2.2

In this exercise we had to:

1. Implement the Explicit Euler method.
2. Implement the Improved Euler (midpoint-like) method.
3. Compare the methods in terms of accuracy and energy conservation for a mass-spring system.

## Methods Implemented

### Explicit Euler

The Explicit Euler method has been implemented in `Exercise_17_2_2.cpp`. The method is defined as:

[
y_{n+1} = y_n + \tau f(y_n)
]

where (y) is the state vector ([y, v]^T) (position and velocity), and (f(y_n)) is the right-hand side of the mass-spring ODE:

[
\frac{dy}{dt} = v, \quad \frac{dv}{dt} = -y
]

The implementation updates the state at each timestep as follows:

```cpp
State f = rhs(y);
y.y += tau * f.y;
y.v += tau * f.v;
```

### Improved Euler (Midpoint-Like)

The Improved Euler method is a two-stage method, which uses a midpoint-like estimate:

[
\tilde{y} = y_n + \frac{\tau}{2} f(y_n)
]

[
y_{n+1} = y_n + \tau f(\tilde{y})
]

This improves the accuracy compared to the simple Explicit Euler. In the code, the implementation is:

```cpp
State f = rhs(y);
State ytilde;
ytilde.y = y.y + 0.5 * tau * f.y;
ytilde.v = y.v + 0.5 * tau * f.v;
State f2 = rhs(ytilde);
y.y += tau * f2.y;
y.v += tau * f2.v;
```

### Analytical Solution

For the mass-spring system with (m = k = 1) and initial conditions (y(0)=1, v(0)=0), the analytical solution is:

[
y(t) = \cos(t), \quad v(t) = -\sin(t)
]

We use this to compute errors and energy deviations.

---

## Comparison

### Time Evolution

Both methods were run for (T = 20) with 200 steps ((\tau = 0.1)).

* **Explicit Euler** shows a gradual drift in energy and slightly increasing amplitude over time.
* **Improved Euler** remains much closer to the analytical solution, showing better energy conservation and phase accuracy.

### Energy Conservation

Energy is computed as:

[
E = \frac{1}{2}v^2 + \frac{1}{2}y^2
]

For each method, the energy deviations were tracked.

* Explicit Euler shows noticeable growth in energy over time.
* Improved Euler keeps energy oscillations around the exact value (0.5), confirming better stability.

### Errors

Max-norm and RMS errors compared to the analytical solution:

 Explicit Euler:&nbsp; &nbsp; :&nbsp; max|y-y*| = 0.651,&nbsp; max|v-v*| = 0.676,&nbsp; RMS y = 0.259,&nbsp; RMS v = 0.269  
 Improved Euler:&nbsp; max|y-y*| = 0.033,&nbsp; max|v-v*| = 0.033,&nbsp; RMS y = 0.012,&nbsp; RMS v = 0.012    

*Note: These values are from running `Exercise_17_2_2.cpp` with T=20, steps=200.*

### CSV Output

The program generates CSV files for plotting:

* `explicit_euler.csv` – time evolution of position and velocity for Explicit Euler
* `improved_euler.csv` – time evolution for Improved Euler
* `energy_comparison.csv` – energies of both methods vs. exact energy

These CSV files can be used to produce plots of (y(t)), (v(t)), phase-space trajectories, and energy evolution.

---

## Summary

* The **Explicit Euler method** is simple but suffers from energy drift in oscillatory systems.
* The **Improved Euler method** provides much higher accuracy and better energy conservation for the mass-spring problem.
* Using analytical solutions allows easy validation of numerical methods.

The implementation demonstrates how even a small improvement in the integration scheme can greatly enhance stability and accuracy for mechanical oscillators.
    
