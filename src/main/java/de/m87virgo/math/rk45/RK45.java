package de.m87virgo.math.rk45;

public class RK45 {

  private final DifferentialEquation equation;
  private final double absTol;
  private final double relTol;
  private final double minStep;
  private final double maxStep;

  public RK45(DifferentialEquation equation, double absTol, double relTol,
              double minStep, double maxStep) {
    this.equation = equation;
    this.absTol = absTol;
    this.relTol = relTol;
    this.minStep = minStep;
    this.maxStep = maxStep;
  }

  /**
   * Integrates y' = f(t,y) from t0 to tEnd with initial condition y0.
   */
  public double integrate(double t0, double y0, double tEnd) {
    double t = t0;
    double y = y0;
    //    double h = (tEnd - t0) / 10.0; // initial step size guess
    double h = minStep; // (tEnd - t0) / 10.0 ist fuer mich sinnbefreit

    int cntIterations = 0;

    while (t < tEnd) {
      if (t + h > tEnd) {
        h = tEnd - t; // last step adjustment
      }
      cntIterations++;
      // Runge-Kutta-Fehlberg coefficients
      double k1 = h * equation.f(t, y);
      double k2 = h * equation.f(t + h / 4.0, y + k1 / 4.0);
      double k3 = h * equation.f(t + 3.0 * h / 8.0, y + 3.0 * k1 / 32.0 + 9.0 * k2 / 32.0);
      double k4 = h * equation.f(t + 12.0 * h / 13.0,
              y + 1932.0 * k1 / 2197.0 - 7200.0 * k2 / 2197.0 + 7296.0 * k3 / 2197.0);
      double k5 = h * equation.f(t + h,
              y + 439.0 * k1 / 216.0 - 8.0 * k2 + 3680.0 * k3 / 513.0 - 845.0 * k4 / 4104.0);
      double k6 = h * equation.f(t + h / 2.0,
              y - 8.0 * k1 / 27.0 + 2.0 * k2 - 3544.0 * k3 / 2565.0
                      + 1859.0 * k4 / 4104.0 - 11.0 * k5 / 40.0);

      // 4th and 5th order estimates
      double y4 = y + 25.0 * k1 / 216.0 + 1408.0 * k3 / 2565.0
              + 2197.0 * k4 / 4104.0 - k5 / 5.0;

      double y5 = y + 16.0 * k1 / 135.0 + 6656.0 * k3 / 12825.0
              + 28561.0 * k4 / 56430.0 - 9.0 * k5 / 50.0 + 2.0 * k6 / 55.0;

      // error estimate
      double err = Math.abs(y5 - y4);
      double tol = absTol + relTol * Math.max(Math.abs(y), Math.abs(y5));
      System.out.println("i:"+ cntIterations+ ", t: " + t + ", y: " + y + ", h: " + h + ", err: " + err + ", tol: " + tol);

      if (err <= tol) {
        // accept step
        t += h;
        y = y5;
      }

      // adaptive step size control
      double safety = 0.9;
      double factor = safety * Math.pow(tol / (err + 1e-15), 0.25);
      factor = Math.max(0.2, Math.min(5.0, factor));

      h *= factor;
      h = Math.max(minStep, Math.min(maxStep, h));
    }
    System.out.println("RK45 iterations: " + cntIterations);
    return y;
  }
}
