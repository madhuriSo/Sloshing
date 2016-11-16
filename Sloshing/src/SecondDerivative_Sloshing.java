
public interface SecondDerivative_Sloshing {
	double getValue(double t, double y, double dy, double[] previousX);
	double getValueByHalfM(double t, double y, double dy, double[] previousX);
}
