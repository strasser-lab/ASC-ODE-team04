# Exercise 18.5 (Pendulum)

In this section, we demonstrate how to implement differential equations without explicitly providing the derivatives (see, for example, Exercise 17.5.1).
To achieve this, we introduce the AutoDiff class, which handles the differentiation automatically.

For the pendulum, the governing second-order ODE is:
$$
\alpha^{\prime \prime} = -\frac{g}{l} \sin (\alpha), \qquad \alpha(t_0) = \pi +0.001, \; \alpha^\prime(t_0) = 0,
$$


Our results are consistent with previous conclusions
- Explicit Euler accumulates numerical errors and introduces artificial energy into the system.
This leads to an overestimation of the amplitude.
- Crank–Nicolson behaves as expected: the oscillations remain stable and physically correct, even with only a tiny perturbation in the initial angle.

## State Space
![](Plots/Pendulum_state_space)
## Pendulum solution
![](Plots/pendulum_Crank_Nicolson.mp4)
![](Plots/pendulum_explicit.mp4)
