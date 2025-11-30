# Exercise 18.4

In this exercise we had to:

1. Add additional useful operators for the AutoDiff class
2. Add some more functions (cos, exp, log, …) for the AutoDiff class.
3. Evaluate and plot Legendre-polynomials up to order 5, in the interval $-1 \le x \le 1$. Evaluate and plot also their derivatives (using AutoDiff).

## Implemented operators/functions (subtask 1, 2)

We implemented some useful operators and functions for the AutoDiff class in the `autodiff.hpp`.

### Arithmetic operators

We have:
$\\ u(x)$ and its derivative $u' \\$
$v(x)$ and its derivative $v' \\$

Then,

#### Addition
Mathematics:

$$
f(x) = u(x) + v(x), \quad f'(x) = u' + v'
$$

Example:

$$
f(x) = x^2 + 3x \Rightarrow f'(x) = 2x + 3
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val + b.val, a.der + b.der);
}
```

#### Subtraction
Mathematics:

$$
f(x) = u(x) - v(x), \quad f'(x) = u' - v'
$$

Example:

$$
f(x) = x^2 - \sin(x) \Rightarrow f'(x) = 2x - \cos(x)
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val - b.val, a.der - b.der);
}
```

#### Multiplication
Mathematics:

$$
f(x) = u(x) \cdot v(x), \quad f'(x) = u'v + uv'
$$

Example:

$$
f(x) = x \cdot \sin(x) \Rightarrow f'(x) = \sin(x) + x\cos(x)
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val * b.val, a.der * b.val + a.val * b.der);
}
```

#### Division
Mathematics:

$$
f(x) = \frac{u(x)}{v(x)}, \quad f'(x) = \frac{u'v - uv'}{v^2}
$$

Example:

$$
f(x) = \frac{x}{x + 1} \Rightarrow f'(x) = \frac{1}{(x + 1)^2}
$$

In the code the implementation is:
```cpp
AutoDiff<T> operator/(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    T bv = b.val;
    return AutoDiff<T>(a.val / bv,
                       (a.der * bv - a.val * b.der) / (bv * bv));
}
```

#### Negative sign
Mathematics:

$$
f(x) = -u(x), \quad f'(x) = -u'
$$

Example:

$$
f(x) = \cos(x) \Rightarrow f'(x) = \sin(x)
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T>& a) {
    return AutoDiff<T>(-a.val, -a.der);
}
```

### Elementary functions
Here we implemented the following useful elementary functions.
We have:
$\\ u(x)$ and its derivative $u' \\$

Then,

#### Sinus
Mathematics:

$$
f(x) = \sin(u(x)) \Rightarrow f'(x) = \cos(u(x)) \cdot u'(x)
$$

Example:

$$
f(x) = \sin(x) \Rightarrow f'(x) = \cos(x)
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> sin(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::sin(a.val), std::cos(a.val) * a.der);
}
```

#### Cosinus
Mathematics:

$$
f(x) = \cos(u(x)) \Rightarrow f'(x) = -\sin(u(x)) \cdot u'(x)
$$

Example:

$$
f(x) = \cos(x) \Rightarrow f'(x) = -\sin(x)
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> cos(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::cos(a.val), -std::sin(a.val) * a.der);
}
```

#### Exponential
Mathematics:

$$
f(x) = e^{u(x)} \Rightarrow f'(x) = e^{u(x)} \cdot u'(x)
$$

Example:

$$
f(x) = e^x \Rightarrow f'(x) = e^x
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> exp(const AutoDiff<T>& a) {
    T ev = std::exp(a.val);
    return AutoDiff<T>(ev, ev * a.der);
}
```

#### Logarithm
Mathematics:

$$
f(x) = \ln(u(x)) \Rightarrow f'(x) = \frac{u'(x)}{u(x)}
$$

Example:

$$
f(x) = \ln(e^x) = x \Rightarrow f'(x) = \frac{1}{e^x} \cdot e^x = 1
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> log(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::log(a.val), a.der / a.val);
}
```

#### Power
Mathematics:

$$
f(x) = (u(x))^c \Rightarrow f'(x) = c \cdot (u(x))^{c-1} \cdot u'(x)
$$

Example:

$$
f(x) = x^3 \Rightarrow f'(x) = 3x^2
$$

In the code the implementation is:
```cpp
template <typename T>
AutoDiff<T> pow(const AutoDiff<T>& a, T exponent) {
    T v = std::pow(a.val, exponent);
    T factor = exponent * std::pow(a.val, exponent - 1);
    return AutoDiff<T>(v, factor * a.der);
}
```

### Scalar operations
Here we implemented the following scalar operations.
We have:
$\\ u(x)$ and its derivative $u' \\$
$c$ as a scalar

Then,

#### Multiplication
Mathematics:

$$
f(x) = c \cdot u(x) \Rightarrow f'(x) = c \cdot u'(x)
$$

$$
f(x) = u(x) \cdot c \Rightarrow f'(x) = u'(x) \cdot c
$$

Example:

$$
f(x) = 5x \Rightarrow f'(x) = 5
$$

In the code the implementation is:
```cpp
template <typename T, typename S>
AutoDiff<T> operator*(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar * a.val, scalar * a.der);
}

template <typename T, typename S>
AutoDiff<T> operator*(const AutoDiff<T>& a, S scalar) {
    return scalar * a;
}
```

#### Addition
Mathematics:

$$
f(x) = c + u(x) \Rightarrow f'(x) = u'(x)
$$

$$
f(x) = u(x) + c \Rightarrow f'(x) = u'(x)
$$

Example:

$$
f(x) = x + 10 \Rightarrow f'(x) = 1
$$

In the code the implementation is:
```cpp
template <typename T, typename S>
AutoDiff<T> operator+(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar + a.val, a.der);
}

template <typename T, typename S>
AutoDiff<T> operator+(const AutoDiff<T>& a, S scalar) {
    return a + scalar;
}
```

#### Subtraction
Mathematics:

$$
f(x) = c - u(x) \Rightarrow f'(x) = -u'(x)
$$

$$
f(x) = u(x) - c \Rightarrow f'(x) = u'(x)
$$

Example:

$$
f(x) = 7 - x \Rightarrow f'(x) = -1
$$

In the code the implementation is:
```cpp
template <typename T, typename S>
AutoDiff<T> operator-(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar - a.val, -a.der);
}

template <typename T, typename S>
AutoDiff<T> operator-(const AutoDiff<T>& a, S scalar) {
    return AutoDiff<T>(a.val - scalar, a.der);
}
```

#### Division
Mathematics:

$$
f(x) = \frac{u(x)}{c} \Rightarrow f'(x) = \frac{u'(x)}{c}
$$

$$
f(x) = \frac{c}{u(x)} \Rightarrow f'(x) = \frac{-c \cdot u'(x)}{u(x)^2}
$$

Example:

$$
f(x) = \frac{x}{2} \Rightarrow f'(x) = \frac{1}{2}
$$

$$
f(x) = \frac{8}{x} \Rightarrow f'(x) = -\frac{8}{x^2}
$$

In the code the implementation is:
```cpp
template <typename T, typename S>
AutoDiff<T> operator/(const AutoDiff<T>& a, S scalar)
{
    return AutoDiff<T>(a.val / scalar, a.der / scalar);
}

template <typename T, typename S>
AutoDiff<T> operator/(S scalar, const AutoDiff<T>& a)
{
    T g  = a.val;
    T g2 = g * g;
    return AutoDiff<T>(scalar / g, -scalar * a.der / g2);
}
```

## Evaluate and plot Legendre-polynomials (and their derivatives) up to order 5, in the interval $-1 \le x \le 1$ (subtask 3)
We implemented this process in `main.cpp`. It evaluates in the interval $-1 \le x \le 1$ the Legendre-polynomials up to order 5 and their derivatives (using AutoDiff). For the plot, it generates a `legendre.csv` file that can be plotted using `plot_legendre.py` program.

### Evaluation
The Legendre polynomials are defined using the following recursive form:

$$
\begin{aligned}
P_{0} = 1 \\
P_{1} = x \\
P_{k}(x) = \frac{(2k - 1)xP_{k-1}(x) - (k - 1)P_{k-2}(x)}{k}, \qquad k \ge 2
\end{aligned}
$$

The program evaluates $P, \, P'$ for 400 ($x$) points in the interval $[-1,1]$. So it computes the values of the Legendre polynomials, and the values of their derivatives for each point.
Then the program writes the results into the `legendre.csv` file in the following format:

$$
x, \quad P_0, \quad P_{0}', \quad P_1, \quad P_{1}', \quad P_2, \quad P_{2}', \quad P_3, \quad P_{3}', \quad P_4, \quad P_{4}', \quad P_5, \quad P_{5}'
$$

For example a line looks like this:

$$
0.5, \quad 1.0, \quad 0.0, \quad 0.5, \quad 1.0, \quad -0.375, \quad 3.0, \quad 0.625, \quad 7.5, \quad -0.625, \quad -15.0, \quad -1.875, \quad -18.75
$$

(So for instance: $P'_{5}(0.5) = -18.75$)

In the end it also prints in the command window, whether the `legendre.csv` was created succesfully.

### Plots
After running the plotting `plot_legendre.py` program, we got the following plot:
![](Plots/legendre_polynomials.png)