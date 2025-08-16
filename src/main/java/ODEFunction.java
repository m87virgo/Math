
/**
 * Scalar ODE right-hand side: y' = f(t, y)
 */
@FunctionalInterface
public interface ODEFunction {
    double f(double t, double y);
}