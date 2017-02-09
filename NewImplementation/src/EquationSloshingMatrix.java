package NewImplementation.src;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class EquationSloshingMatrix {
    //this parameter is not configurable
	public static double g=9.81;

    //Parameters of liquid
	public static double RHO=1000.0;
	public static double MU=0.00152;
	public static double h=0.05;
	public static double ZETA=0.07;

    //Parameters of liquid with respect to time
	public static int T_MAX=60;
	public static double DELTA_T=1.0/256;

    //Parameter of the tank
	public static double L=1.0;
	public static int N=50;
	double l=(L/N);

    //Parameters for external vibrations
	public static double A=0.0035;
	public static double F_REQ=0.386;

    double m = 0; 		// to be initialized
    double cc = 0 ;		//damping coeff
    double kb = 0;		// base spring
    double cb = 0;		// connecting damper
    double OMEGA = 0;

    double 	ETA=0;

    public EquationSloshingMatrix()
    {
        this.cc = ZETA / l ;
        this.OMEGA = 2 * Math.PI * F_REQ;
        this.cb = MU * l * Math.sqrt((OMEGA * RHO)/(2*MU));
        this.kb = - OMEGA * cb;
        this.m = RHO * l * h;
    }

    // Calculate alpha + beta
    public double getAlphaPlusBeta (double eta,double x0_dot, double x1_dot, double x0,double x1) {
    	double alpha, beta;
        alpha=0.5*  RHO * g *((eta*eta)+(2*eta*h));
        beta=((1.0/3) *RHO*(eta+h)*(eta+h)*(2*h*l*(x1_dot-x0_dot)*(x1_dot-x0_dot))/((x1-x0+l)*(x1-x0+l)*(x1-x0+l)));
        return alpha + beta;
    }
    
    public double getAlphaPlusBetaTemp (double eta,double x0_dot, double x1_dot, double x0,double x1) {
    	double alpha, beta;
        alpha=0.5*  RHO * g *((eta*eta)+(2*eta*h));
        beta=((1.0/3) *RHO*(eta+h)*(eta+h)*(2*h*l*(x1_dot-x0_dot)*(x1_dot-x0_dot))/((x1-x0+l)*(x1-x0+l)*(x1-x0+l)));
        return alpha + beta;
    }
    
    // Calculate ETA
    public double getEta(double x1, double x0){
        return h*((l/(x1-x0+l))-1);
    }

    // Calculate Fe
    public double getFe(double t){
        return -m* (OMEGA*OMEGA)* A *(Math.sin(OMEGA*t));
    }
    

    //get value of F(x,xdot)
    //               0      1         2              3        4          5            6       7          8
    //preResult[]= {x0,   x0_dot,    x0_doubleDot,  x1,    x1_dot,   x1_doubleDot,   x2,   x2_dot,   x2_doubleDot}

	
	public double getValue(double t, double[] previousX) {
		double f = 0 ,alphaPlusBeta1 = 0, alphaPlusBeta2=0,fe;
        double x0=previousX[0],x0_dot=previousX[1];
        double x1=previousX[3],x1_dot=previousX[4];
        double x2=previousX[6],x2_dot=previousX[7];
        ETA=getEta(x1,x0);
        alphaPlusBeta1=getAlphaPlusBeta(ETA,x0_dot,x1_dot,x0,x1);
        ETA=getEta(x2,x1);
        alphaPlusBeta2=getAlphaPlusBeta(ETA,x1_dot,x2_dot,x1,x2);
        fe=getFe(t);
        f=(alphaPlusBeta1 - alphaPlusBeta2) + cc * x0_dot - (((2 * cc) + cb) * x1_dot) + (cc * x2_dot) - (kb * x1) + fe;
        return f;
	}
	
	// Calculate Gamma i
	public double getGammaCurrent (double[] previousX) {
		double x0=previousX[0]; // x_i-1
        double x1=previousX[3]; // x_i
        ETA=getEta(x1,x0);
        return (1.0/3)*RHO*h*l*((ETA+h)*(ETA+h))/((x1-x0+l)*(x1-x0+l));
	}
	
	// Calculate Gamma i+1
		public double getGammaNext (double[] previousX) {
	        double x1=previousX[3]; // x_i
	        double x2=previousX[6]; // x_i+1
	        ETA=getEta(x2,x1);
	        double gamma = (1.0/3)*RHO*h*l*((ETA+h)*(ETA+h))/((x2-x1+l)*(x2-x1+l));
	        return gamma;
		}
		
	//Load values from configuration file
		
	public  void parseFile(String filename) 
	{
//			public static void main(String[] args) {
			Properties prop = new Properties();
			BufferedReader input = null;
			try 
			{
				input = new BufferedReader(new FileReader(filename));		
				// load a properties file
				prop.load(input);
				
				 for ( String line; null != ( line = input.readLine() );  ) 
			        {
			        	if ("".equals( line.trim() ) ||line.startsWith( "*" ))
			                return;
			        }

				
				// get the property value and print it out
//				System.out.println(prop.getProperty("g"));
//				System.out.println(prop.getProperty("RHO"));
//				System.out.println(prop.getProperty("MU"));
//				System.out.println(prop.getProperty("h"));
//				System.out.println(prop.getProperty("ZETA"));
//				System.out.println(prop.getProperty("T_MAX"));
//				System.out.println(prop.getProperty("DELTA_T"));
//				System.out.println(prop.getProperty("L"));
//				System.out.println(prop.getProperty("N"));
//				System.out.println(prop.getProperty("A"));
//				System.out.println(prop.getProperty("F_REQ"));
				g=Double.parseDouble(prop.getProperty("g"));
				RHO=Double.parseDouble(prop.getProperty("RHO"));
				MU=Double.parseDouble(prop.getProperty("MU"));
				h=Double.parseDouble(prop.getProperty("h"));
				ZETA=Double.parseDouble(prop.getProperty("ZETA"));
				T_MAX=Integer.parseInt(prop.getProperty("T_MAX"));
				DELTA_T=Double.parseDouble(prop.getProperty("DELTA_T"));
				L=Double.parseDouble(prop.getProperty("L"));
				N=Integer.parseInt(prop.getProperty("N"));
				A=Double.parseDouble(prop.getProperty("A"));
				F_REQ=Double.parseDouble(prop.getProperty("F_REQ"));
				//System.out.println(DELTA_T);
			}  
			catch (FileNotFoundException e) 
			{
				e.printStackTrace();
			} 
			catch (IOException ex) 
			{
				ex.printStackTrace();
			} 
			finally 
			{
				if (input != null)
				{
					try 
					{
						input.close();
					} 
					catch (IOException e) 
					{
						e.printStackTrace();
					}
				}
			}
	}
}
