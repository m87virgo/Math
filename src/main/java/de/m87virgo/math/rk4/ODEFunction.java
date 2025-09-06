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
    final class RK4 {

        private RK4() { }

        /**
         * Advance one RK4 step.
         * @param f   right-hand side f(t, y)
         * @param t   current time
         * @param y   current state
         * @param h   step size
         * @return    y_{n+1}
         */
        public static double step(ODEFunction f, double t, double y, double h) {
            double k1 = f.f(t, y);
            double k2 = f.f(t + 0.5 * h, y + 0.5 * h * k1);
            double k3 = f.f(t + 0.5 * h, y + 0.5 * h * k2);
            double k4 = f.f(t + h,       y + h * k3);
            return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

        /**
         * Integrate from t0 to tEnd in nSteps uniform steps.
         * @return array of y values of length nSteps+1 (including y0)
         */
        public static double[] solve(ODEFunction f, double t0, double y0, double tEnd, int nSteps) {
            if (nSteps <= 0) throw new IllegalArgumentException("nSteps must be > 0");
            double[] ys = new double[nSteps + 1];
            ys[0] = y0;
            double t = t0;
            double y = y0;
            double h = (tEnd - t0) / nSteps;
            for (int i = 1; i <= nSteps; i++) {
                y = step(f, t, y, h);
                t += h;
                ys[i] = y;
            }
            return ys;
        }
    }
}