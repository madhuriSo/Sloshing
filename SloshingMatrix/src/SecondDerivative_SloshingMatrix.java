
public interface SecondDerivative_SloshingMatrix {
	// get value of F(x,xdot)
	double getValue(double t, double y, double dy, double[] previousX);
}
