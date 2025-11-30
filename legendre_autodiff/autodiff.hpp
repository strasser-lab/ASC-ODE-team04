// autodiff.hpp
#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <array>
#include <cmath>
#include <iostream>

// ================================================================
// Simple forward-mode AutoDiff for ONE variable (perfect for Legendre demo)
// ================================================================

template <typename T = double>
class AutoDiff {
public:
    T val;      // function value
    T der;      // derivative (only one variable → single number)

    // constructors
    AutoDiff() : val(0), der(0) {}
    AutoDiff(T v) : val(v), der(0) {}
    AutoDiff(T v, T d) : val(v), der(d) {}

    // create the independent variable x (seed derivative = 1)
    static AutoDiff variable(T v) { return AutoDiff(v, T(1)); }
};

// --------------------- Arithmetic operators ---------------------

// addition
template <typename T>
AutoDiff<T> operator+(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val + b.val, a.der + b.der);
}

// subtraction
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val - b.val, a.der - b.der);
}

// multiplication
template <typename T>
AutoDiff<T> operator*(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    return AutoDiff<T>(a.val * b.val, a.der * b.val + a.val * b.der);
}

// division
template <typename T>
AutoDiff<T> operator/(const AutoDiff<T>& a, const AutoDiff<T>& b) {
    T bv = b.val;
    return AutoDiff<T>(a.val / bv,
                       (a.der * bv - a.val * b.der) / (bv * bv));
}

// negative sign
template <typename T>
AutoDiff<T> operator-(const AutoDiff<T>& a) {
    return AutoDiff<T>(-a.val, -a.der);
}

// --------------------- Elementary functions ---------------------

// sinus
template <typename T>
AutoDiff<T> sin(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::sin(a.val), std::cos(a.val) * a.der);
}

// cosinus
template <typename T>
AutoDiff<T> cos(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::cos(a.val), -std::sin(a.val) * a.der);
}

// exponential
template <typename T>
AutoDiff<T> exp(const AutoDiff<T>& a) {
    T ev = std::exp(a.val);
    return AutoDiff<T>(ev, ev * a.der);
}

// logarithm
template <typename T>
AutoDiff<T> log(const AutoDiff<T>& a) {
    return AutoDiff<T>(std::log(a.val), a.der / a.val);
}

// power
template <typename T>
AutoDiff<T> pow(const AutoDiff<T>& a, T exponent) {
    T v = std::pow(a.val, exponent);
    T factor = exponent * std::pow(a.val, exponent - 1);
    return AutoDiff<T>(v, factor * a.der);
}

// --------------------- Scalar operations ---------------------
// also allow scalar operations (very convenient)

// multiplication
//    f(x) = c · g(x)    →    f' = c · g'
template <typename T, typename S>
AutoDiff<T> operator*(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar * a.val, scalar * a.der);
}

template <typename T, typename S>
AutoDiff<T> operator*(const AutoDiff<T>& a, S scalar) {
    return scalar * a;
}

// addition
//    f(x) = c + g(x)    →    f' = g'
//    f(x) = g(x) + c    →    f' = g'
template <typename T, typename S>
AutoDiff<T> operator+(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar + a.val, a.der);
}

template <typename T, typename S>
AutoDiff<T> operator+(const AutoDiff<T>& a, S scalar) {
    return a + scalar;
}

// subtraction
//    f(x) = c - g(x)    →    f' = -g'
//    f(x) = g(x) - c    →    f' = g'
template <typename T, typename S>
AutoDiff<T> operator-(S scalar, const AutoDiff<T>& a) {
    return AutoDiff<T>(scalar - a.val, -a.der);
}

template <typename T, typename S>
AutoDiff<T> operator-(const AutoDiff<T>& a, S scalar) {
    return AutoDiff<T>(a.val - scalar, a.der);
}

// division
// f(x) = g(x) / c     →     f' = g'/c
// f(x) = c / g(x)     →     f' = -c * g' / g²
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

#endif // AUTODIFF_HPP