package de.m87virgo.math.rk45;

@FunctionalInterface
public interface DifferentialEquation {
  /**
   * Computes the derivative dy/dt = f(t, y).
   *
   * @param t current time
   * @param y current state
   * @return derivative f(t, y)
   */
  double f(double t, double y);
}
