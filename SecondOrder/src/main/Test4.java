package main;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public class Test4 {
    public static void main(String args[]) {
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

        System.out.println("TEST4, equation of motion");
        System.out.println(" t,s             y,m             v,m/s               a,m/s2");
        System.out.println("--------------------------------------------------------------");

        System.out.println("\n Time0    \t Position0    \t Velocity   \t   acceleration");
        System.out.println("\t" +integratorRKN.t0+ " \t   "+integratorRKN.y0+"      \t   " +integratorRKN.dy0+"       \t  " +integratorRKN.d2y0);

        double t=t0, y = y0, v=v0, a=a0;

        while (t < tmax) {
            double[] newValues=integratorRKN.Step(t, y, v, a);
            System.out.print("\t"+newValues[0]+"\t"+newValues[1]+"\t"+newValues[2]+"\t"+newValues[3]);
            t=newValues[0];
            y=newValues[1];
            v=newValues[2];
            a=newValues[3];
            System.out.print("\n");

        }
    }
}
