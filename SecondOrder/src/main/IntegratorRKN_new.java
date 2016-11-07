
public class IntegratorRKN_new {
	SecondDerivative secondDerivative;
    double t0, y0, dy0, h;
    double d2y0;

    public IntegratorRKN_new (SecondDerivative secondDerivative, double t0, 
        double y0, double dy0, double h, double d2y0)
    {
       this.secondDerivative = secondDerivative;
       this.t0 = t0;
       this.y0 = y0;
       this.dy0 = dy0;
       this.h = h;
       this.d2y0 = d2y0;
    }

    // Different from previous code, t0 will not be updated in this method.
    public IntegratorRKN_new step(double t, double y, double dy,  double d2y)
    {
       double h2 = h * h; 

       double k1 = d2y0;
       double k2 = secondDerivative.getValue(
           t0 + h/2, y0 + h/2 * dy0 + h2/8 * k1, dy0 + h/2 * k1);
       double k3 = secondDerivative.getValue(
           t0 + h/2, y0 + h/2 * dy0 + h2/8 * k2, dy0 + h/2 * k2);
       double k4 = secondDerivative.getValue(
           t0 + h,   y0 +   h * dy0 + h2/2 * k3, dy0 +   h * k3);

       //t = t0 + h;
       y = y0 + h * dy0 + h2/6 * (k1 + k2 + k3);
       dy = dy0 + h/6 * (k1 + 2 * k2 + 2 * k3 + k4);
       d2y = secondDerivative.getValue(t, y, dy);

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
