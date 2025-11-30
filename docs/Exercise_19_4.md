**Problem Description**
This exercise aims to compare different Runge-Kutta (RK) methods in terms of numerical accuracy and convergence when solving a first-order ordinary differential equation (ODE).  
The ODE considered is:  

$$
\begin{align}
    y′(t)=−y(t), \\
    y(0)=1  
\end{align}
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
\begin{align}
    k_1​=f(t_n​,y_n​), \\
    k_2​=f(t_n​ +h/2,y_n​+h k_1​/2), \\ 
    y_{n+1}​=y_n+h k_2​  
\end{align}    
$$  

Characteristics: second-order convergence, explicit computation.  
RK4  
+ Type: Explicit fourth-order Runge-Kutta  
+ Step formula:  
  
$$  
\begin{align}
    k_1​ = f(t_n​,y_n​), \\
    k_2​ = f(t_n​+h/2,y_n​+h k_1​ /2), \\
    k_3​ = f(t_n​+h/2,y_n​+h k_2​ /2), \\
    k_4​=f(t_n​ + h,y_n​+h k_3​), \\
    y_{n+1}​=y_n​+6/h​(k_1​+2 k_2​+2 k_3​+k_4) \\
\end{aligned}
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
