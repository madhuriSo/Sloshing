
public class FlyingBall_2DArray implements SecondDerivative_2DArray {
	double dB;    // Ball diameter, m
    double roB;   // Ball density, kg/m^3
    double cD;    // Drag coefficient
    double g;     // Gravitational acceleration, m/s^2
    double roA;   // Air density
    double mB;    // Ball mass, kg
    double aB;    // Ball cross-section area, m^2

    public FlyingBall_2DArray()
    {
       dB = 0.1;
       roB = 600;
       cD = 0.1;
       g = 9.81;
       roA = 1.29;

       double v = Math.PI * Math.pow(dB, 3) / 6;
       mB = roB * v;
       aB = 0.25 * Math.PI * dB * dB;
    }

    public double getValue(double t, double y, double v, double[] xValues)
    {
       double f = -mB * g -Math.signum(v) * 0.5 * cD * roA * v * v * aB;
       double d2y = f / mB;
       return d2y;
    }
}
