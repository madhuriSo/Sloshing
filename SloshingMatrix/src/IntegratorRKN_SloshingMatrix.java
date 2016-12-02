
public class IntegratorRKN_SloshingMatrix {
	SecondDerivative_SloshingMatrix secondDerivative;
    double t0, y0, dy0, h;
    double d2y0;
    double [] previousX;

    public IntegratorRKN_SloshingMatrix (SecondDerivative_SloshingMatrix secondDerivative, double t0, 
        double y0, double dy0, double h, double d2y0)
    {
       this.secondDerivative = secondDerivative;
       this.t0 = t0;
       this.y0 = y0;
       this.dy0 = dy0;
       this.h = h;
       this.d2y0 = d2y0;
       previousX = new double[9]; // to store xi-1, xi, xi+1 etc., total 9 values
    }

    // Different from previous code, t0 will not be updated in this method.
    public IntegratorRKN_SloshingMatrix step(double t, double y, double dy,  double d2y, double[] previousX)
    {
       double h2 = h * h; 

       double k1 = d2y0;
       double k2 = secondDerivative.getValue(
           t0 + h/2, y0 + h/2 * dy0 + h2/8 * k1, dy0 + h/2 * k1, previousX);
       double k3 = secondDerivative.getValue(
           t0 + h/2, y0 + h/2 * dy0 + h2/8 * k2, dy0 + h/2 * k2, previousX);
       double k4 = secondDerivative.getValue(
           t0 + h,   y0 +   h * dy0 + h2/2 * k3, dy0 +   h * k3, previousX);

       //t = t0 + h;
       y = y0 + h * dy0 + h2/6 * (k1 + k2 + k3);
       dy = dy0 + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
       d2y = secondDerivative.getValue(t, y, dy, previousX);

       //t0 = t;
       y0 = y;
       dy0 = dy;
       d2y0 = d2y;

       return this;
}
    
    public void updateT0 () {
    	t0 += h;
    }
    
}
