package main;
import java.lang.*;

/**
 * Created by Madhuri on 11/4/16.
 */

public class FlyingBall implements SecondDerivative {
    double dB;    // Ball diameter, m
    double roB;   // Ball density, kg/m^3
    double cD;    // Drag coefficient
    double g;     // Gravitational acceleration, m/s^2
    double roA;   // Air density
    double mB;    // Ball mass, kg
    double aB;    // Ball cross-section area, m^2

    public FlyingBall()
    {
        dB = 0.1;
        roB = 600;
        cD = 0.1;
        g = 9.81;
        roA = 1.29;


        double v = (Math.PI * Math.pow(dB,3))/6;
        mB = roB * v;
        aB = 0.25 * Math.PI * dB * dB;
    }


    @Override
    public double GetValue(double t, double y, double v) {
        double f = (-(mB * g) -(Math.signum(v)* 0.5 * cD * roA * v * v * aB));
        double d2y = f / mB;
        return d2y;

    }
}
