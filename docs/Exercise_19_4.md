**Problem Description**
This exercise aims to compare different Runge-Kutta (RK) methods in terms of numerical accuracy and convergence when solving a first-order ordinary differential equation (ODE).  
The ODE considered is:  

$$
\begin{aligned}
    y′(t)=−y(t), \\
    y(0)=1  
\end{aligned}
$$  

Its analytical solution is:  

$$  
    y(t)=e^{−t}
$$  

The numerical methods to be compared are:  
  1. Explicit RK2 (Midpoint)  
  2. Explicit RK4  
  3. Implicit Gauss-Legendre 2-stage (s=2)  
  4. Implicit Radau IIA 2-stage (s=2)

We compute the approximate values and errors at each time step  
  
**Principles of the Four Methods**  
RK2 (Midpoint)  
+ Type: Explicit second-order Runge-Kutta  
+ Step formula:  
  
$$
\begin{aligned}
    k_1​=f(t_n​,y_n​), \\
    k_2​=f(t_n​ +h/2,y_n​+h k_1​/2), \\ 
    y_{n+1}​=y_n+h k_2​  
\end{aligned}
$$  

Characteristics: second-order convergence, explicit computation.  
RK4  
+ Type: Explicit fourth-order Runge-Kutta  
+ Step formula:  
  
$$
$$  k_1​=f(t_n​,y_n​), $$
$$  k_2​=f(t_n​+h/2,y_n​+h k_1​/2), $$
$$  k_3​=f(t_n​+h/2,y_n​+h k_2​/2), $$
$$  k_4​=f(t_n​ + h,y_n​+h k_3​), $$
$$  y_{n+1}​=y_n​+6/h​(k_1​+2 k_2​+2 k_3​+k_4) $$
$$  

Caracteristics: fourth-order convergence, explicit, higher accuracy than RK2.  
Gauss-Legendre 2-stage  
+ Type: Implicit two-stage Gauss-Legendre Runge-Kutta  
+ Uses two nodes (s=2) in a high-order implicit integration formula  
+ Each step requires solving a nonlinear system (fixed-point iteration or Newton iteration)  
+ Characteristics: A-stable, second-order accuracy, implicit computation.

Radau IIA 2-stage  
+ Type: Implicit Radau IIA Runge-Kutta  
+ Two nodes, implicit method  
+ Suitable for stiff ODEs, A-stable and second-order accurate  
+ Step requires solving a nonlinear system.

**Error Comparison and Results**  
Calculation conditions:  
+ Final time 𝑇 = 1.0  
+ Number of steps 𝑁 = 10 (step size ℎ=0.1)

Sample numerical results:  
| t        | RK2        | error      | RK4        | error      | GL2        | error      | Radau2     | error      |
|----------|------------|------------|------------|------------|------------|------------|------------|------------|
| 0.0      | 1.00000000 | 0.00000000 | 1.00000000 | 0.00000000 | 1.00000000 | 0.00000000 | 1.00000000 | 0.00000000 |
| 0.1      | 0.90500000 | 0.00016258 | 0.90483750 | 0.00000008 | 0.90483743 | 0.00000001 | 0.90483619 | 0.00000122 |
| 0.2      | 0.81902500 | 0.00029425 | 0.81873090 | 0.00000015 | 0.81873078 | 0.00000002 | 0.81872854 | 0.00000222 |
| 0.3      | 0.74121762 | 0.00039940 | 0.74081842 | 0.00000020 | 0.74081825 | 0.00000003 | 0.74081521 | 0.00000301 |
| 0.4      | 0.67080195 | 0.00048190 | 0.67032029 | 0.00000024 | 0.67032008 | 0.00000004 | 0.67031642 | 0.00000363 |
| 0.5      | 0.60707577 | 0.00054511 | 0.60653093 | 0.00000027 | 0.60653070 | 0.00000004 | 0.60652656 | 0.00000410 |
| 0.6      | 0.54940357 | 0.00059193 | 0.54881193 | 0.00000030 | 0.54881168 | 0.00000005 | 0.54880718 | 0.00000446 |
| 0.7      | 0.49721023 | 0.00062492 | 0.49658562 | 0.00000031 | 0.49658535 | 0.00000005 | 0.49658060 | 0.00000470 |
| 0.8      | 0.44997526 | 0.00064629 | 0.44932929 | 0.00000033 | 0.44932901 | 0.00000005 | 0.44932410 | 0.00000486 |
| 0.9      | 0.40722761 | 0.00065795 | 0.40656999 | 0.00000033 | 0.40656971 | 0.00000005 | 0.40656471 | 0.00000495 |
| 1.0      | 0.36854098 | 0.00066154 | 0.36787977 | 0.00000033 | 0.36787949 | 0.00000005 | 0.36787446 | 0.00000498 |

Observations:  
  1. RK2 exhibits larger errors, deviating noticeably from the exact solution.  
  2. RK4 has very small errors, demonstrating clear advantage of fourth-order convergence.  
  3. Gauss-Legendre 2-stage and Radau IIA 2-stage have very small errors, comparable to RK4, and are A-stable implicit methods.  
  4. Reducing step size decreases the errors of RK4, GL2, and Radau2 rapidly, while RK2 decreases more slowly.  
  
**Conclusion**  
a.Explicit methods:  
+ RK2 has limited accuracy, suitable for non-stiff ODEs with low precision requirements.  
+ RK4 is fourth-order accurate, computationally efficient, and suitable for general non-stiff problems.
  
b.Implicit methods:  
+ Gauss-Legendre and Radau IIA are more stable for stiff problems (A-stable).  
+ Accuracy is comparable to RK4, and larger step sizes can be used without divergence.
  
c.Summary:  
+ For non-stiff ODEs with high accuracy requirements, RK4 is the simplest and most effective choice.  
+ For stiff ODEs or when larger step sizes are desired, implicit methods (Gauss-Legendre / Radau IIA) are more appropriate.  
