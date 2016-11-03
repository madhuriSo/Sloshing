package main;


public class RangeKutta {

    // The number of steps to use in the interval
    public static final double N = 2;


    // The derivative dx/dt at a given value of t and x.
    public static double deriv(double t, double x)
    {
        double a,b,c;
        a=Math.exp(-t);
        b=0.4*x;
        c=((3*(Math.exp(-t)))-(0.4*x));
        return (3*(Math.exp(-t))-(0.4*x));
    }


    // The `main' method does the actual computations
    public static void main(String[] argv)
    {
        // `h' is the size of each step.
        double h = 1.5;
        double A, B, C, D;
        double t, x;
        int i;


        // Computation by 4th order Runge-Kutta

        // Initialize x
        x = 5;

        for (i=0; i<N; i++)
        {
            // Step through, updating t
            t = i * h;

            // Computing all of the trial values
            A = h * deriv(t, x);
            B = h * deriv(t + h/2, x + A/2);
            C = h * deriv(t + h/2, x + B/2);
            D = h * deriv(t + h, x + C);

            // Incrementing x
           x += A/6 + B/3+ C/3 + D/6;
            System.out.println(" Value of X: "+x+" at t= "+i);


        }

        // Print out the result that we get.

        System.out.println("Using 4th order Runge-Kutta "
                + "The value at t=3 is:");
        System.out.println(x);


    }
}

