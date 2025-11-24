package de.m87virgo.math.rk4;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

/**
 * Unit tests comparing numerical RK4 with analytical solutions.
 */
public class RK4Test {

    @Test
    void testExpGrowth_yPrimeEqualsY() {
        // y' = y, y(0) = 1 -> y(t) = e^t
        ODEFunction f = (t, y) -> y;

      double t0 = 0.0;
        double y0 = 1.0;
        double tEnd = 1.0;     // compare at t = 1

        int nSteps = 100;      // sufficiently small h
        double[] ys = RK4.solve(f, t0, y0, tEnd, nSteps);
        double yNum = ys[ys.length - 1];
        double yExact = Math.exp(1.0);

        double absErr = Math.abs(yNum - yExact);
        assertTrue(absErr < 1e-8, "RK4 solution should closely match e at t=1; error=" + absErr);
    }

    @Test
    void testGlobalErrorOrderApproximatelyFourth() {
        // Same ODE, compare error reduction when halving h
        ODEFunction f = (t, y) -> y;
        double t0 = 0.0;
        double y0 = 1.0;
        double tEnd = 1.0;

        int n1 = 25;            // coarse
        int n2 = 50;            // half step size
        double exact = Math.exp(1.0);

        double[] y1 = RK4.solve(f, t0, y0, tEnd, n1);
        double[] y2 = RK4.solve(f, t0, y0, tEnd, n2);
        double e1 = Math.abs(y1[y1.length - 1] - exact);
        double e2 = Math.abs(y2[y2.length - 1] - exact);

        // For RK4, halving h should reduce error by about 2^4 = 16.
        // Allow some slack because constant factors differ.
        assertTrue(e1 > 10 * e2, "Error should drop ~16x; got e1/e2 = " + (e1 / e2));
    }

    @Test
    void testLinearODE_yPrimeEqualsMinus2yPlusT() {
        // y' = -2y + t, y(0) = 1
        // Analytical solution: y(t) = 1.25 * e^{-2t} + 0.5 t - 0.25
        ODEFunction f = (t, y) -> -2.0 * y + t;
        double t0 = 0.0;
        double y0 = 1.0;
        double tEnd = 2.0;
        int nSteps = (int)1e5;
        double[] ys = RK4.solve(f, t0, y0, tEnd, nSteps);
        double yNum = ys[ys.length - 1];
        double yExact = 1.25 * Math.exp(-2.0 * tEnd) + 0.5 * tEnd - 0.25;
        assertEquals(yExact, yNum, 1e-7);
    }
}