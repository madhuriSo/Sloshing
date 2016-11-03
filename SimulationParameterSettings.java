

/**
 * @author Richa Gupta
 *Please refer to slosh-4.pptx slide 17
 */

public class SimulationParameterSettings 
{
	//this parameter is not configurable
	public static final double g=9.81;
	
	//Parameters of liquid
	public static final int RHO=1000;
	public static final double MU=0.00152;
	public static final double h=0.05;
	public static final double ZETA=0.07;
	
	//Parameters of liquid with respect to time
	public static final int T_MAX=8;
	public static final double DELTA_T=1/1024;
	
	//Parameter of the tank
	public static final double L=1.0;
	public static final double N=100;
	public static final double l=(L/N);
	
	//Parameters for external vibrations
	public static final double A=0.0035;
	public static final double F_REQ=0.386;
	
	
}
