package de.m87virgo.math.rk45;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

public class RK45Test {

  @Test
  public void testExponentialGrowth() {
    // DGL: y' = y, Lösung: y(t) = exp(t), Startwert: y(0) = 1
    DifferentialEquation eq = (t, y) -> y;
    RK45 solver = new RK45(eq, 1e-8, 1e-8, 1e-6, 0.5);

    double yEnd = solver.integrate(0.0, 1.0, 1.0);
    double exact = Math.exp(1.0);

    assertEquals(exact, yEnd, 1e-6, "RK45 should approximate exp(1)");
  }

  @Test
  public void testSinCosRelation() {
    // DGL: y' = cos(t), Lösung: y(t) = sin(t), Startwert: y(0) = 0
    DifferentialEquation eq = (t, y) -> Math.cos(t);
    RK45 solver = new RK45(eq, 1e-8, 1e-8, 1e-6, 0.5);

    double yEnd = solver.integrate(0.0, 0.0, Math.PI / 2);
    double exact = Math.sin(Math.PI / 2);

    assertEquals(exact, yEnd, 1e-6, "RK45 should approximate sin(pi/2) = 1");
  }

  /**
   * Test mit einer funktion die einen Wendepunkt hat.
   *  DGL: y' = t, Lösung: y(t) = 0.5 * t^2, Startwert: y(0) = 0
   *  Dieser Test prüft, ob der RK45-Integrator auch bei einem Wendepunkt
   *  (hier bei ( t=0 )) die Schrittweite korrekt anpasst und das Ergebnis präzise bleibt.
   */
  @Test
  public void testParabolaWithInflectionPoint() {
    DifferentialEquation eq = (t, y) -> t;
    RK45 solver = new RK45(eq, 1e-8, 1e-8, 1e-8, 0.1);

    double t0 = -2.0;
    double tEnd = 2.0;
    double y0 = 0.0;
    double yEnd = solver.integrate(t0, y0, tEnd);
    double exact = 0.5 * tEnd * tEnd;

    assertEquals(exact, yEnd, 1e-6, "RK45 sollte Wendepunkt bei t=0 korrekt behandeln");
  }
}
