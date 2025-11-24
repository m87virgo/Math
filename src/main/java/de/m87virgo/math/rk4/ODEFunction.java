package de.m87virgo.math.rk4;

/**
 * Scalar ODE right-hand side: y' = f(t, y)
 */
@FunctionalInterface
public interface ODEFunction {
    double f(double t, double y);

    /**
     * Classic 4th-order Runge–Kutta for scalar ODEs.
     */

}