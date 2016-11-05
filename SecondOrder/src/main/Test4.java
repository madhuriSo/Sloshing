package main;


import java.io.Console;

/**
 * Created by Madhuri on 11/4/16.
 */
public class Test4 {
    FlyingBall ball = new FlyingBall();
    SecondDerivative acceleration = (SecondDerivative) ball;

    double t0 = 0;
    double tmax = 10;
    double y0 = 0;
    double h = 0.5;
    double v0 = 50;
    double a0 = acceleration.GetValue(t0, y0, v0);
    IntegratorRKN integratorRKN = new IntegratorRKN(acceleration, t0, y0, v0,
            h, a0);

    System.out.println("Name entered :");
    Console.WriteLine("TEST4, equation of motion");
    Console.WriteLine(" t,s     y,m   v,m/s  a,m/s2");
    Console.WriteLine("----------------------------");
    Console.WriteLine("{0,4:F1}{1,8:F2}{2,8:F2}{3,8:F2}", t0, y0, v0, a0);

    double t, y, v, a;
    while(t < tmax)
    {
        integratorRKN.Step(t, y, v, a);
        Console.WriteLine("{0,4:F1}{1,8:F2}{2,8:F2}{3,8:F2}", t, y, v, a);
    }
    Console.WriteLine();
}
