import java.io.FileWriter;
import java.io.IOException;

public class SolveEquation_Method1 {
	public static void main (String[] args) {
		double t = 0;
		double deltaT = 1.0 / 256;
		double T_MAX = 0.05;
		int N = 20;
	
		double[] xPrevious = new double[N+1];
		double[] xDotPrevious = new double[N+1];
		double[] xDDotPrevious = new double[N+1];
		
		double[] x = new double[N+1]; // current
		double[] xDot = new double[N+1]; // current
		double[] xDDot = new double[N+1]; // current
	
		double[] xPrime = new double[N+1];
		double[] xDotPrime = new double[N+1];
		double[] xDDotPrime = new double[N+1];
	
		double[] A_x = new double[N+1];
		double[] A_xDot = new double[N+1];
		double[] B_x = new double[N+1];
		double[] B_xDot = new double[N+1];
		double[] C_x = new double[N+1];
		double[] C_xDot = new double[N+1];
		double[] D_x = new double[N+1];
		double[] D_xDot = new double[N+1];
		
		double[] previousX = new double[9];
		double[] previousXPrime = new double[9];
		
		double[][] allXResult = new double[(int)(T_MAX/deltaT)+1][3*(N+1)+1];
		double[][] etas = new double[(int)(T_MAX/deltaT)+1][N+1];
		double[][] allTs = new double[(int)(T_MAX/deltaT)+1][1];
		int rowCounter = 0;
		
		EquationSloshing sloshingEquation = new EquationSloshing();
	
		while (t <= T_MAX) {
			allTs[rowCounter][0] = t;
			System.out.println(t);
			
			//get A_x
			for (int i = 0; i < N + 1; i++) {
				A_x[i] = xDotPrevious[i];
			}
			
			//compute A_xDot
			A_xDot[0] = 0;
			for (int j = 1; j < N; j++) {
				// get previousX
				previousX[0] = xPrevious[j-1];
				previousX[1] = xDotPrevious[j-1];
				previousX[2] = xDDotPrevious[j-1];
				previousX[3] = xPrevious[j];
				previousX[4] = xDotPrevious[j];
				previousX[5] = xDDotPrevious[j];
				previousX[6] = xPrevious[j+1];
				previousX[7] = xDotPrevious[j+1];
				previousX[8] = xDDotPrevious[j+1];
				
				//System.out.println("j=" + j + "  "+ previousX[0] + "  " + previousX[1] + "  " + previousX[2] + "  " +
				//                                            previousX[3] + "  " + previousX[4] + "  " + previousX[5] + "  " +
				//		                                    previousX[6] + "  " + previousX[7] + "  " + previousX[8]);
				
				A_xDot[j] = sloshingEquation.getValue(t, previousX);
				//System.out.println("j=" + j + "  " + "A_xDot=" + A_xDot[j]);
			}
			A_xDot[N] = 0;
			
			// compute B_x, B_xDot
			// compute xPrime, xDotPrime
			for (int k = 0; k < N; k++) {
				xPrime[k] = xPrevious[k] + deltaT / 2 * A_x[k];
				xDotPrime[k] = xDotPrevious[k] + deltaT / 2 * A_xDot[k];
			}
			// get B_x
			for (int k = 0; k < N; k++) {
				B_x[k] = xDotPrime[k];
			}
			// compute B_xDot
			B_xDot[0] = 0;
			for (int j = 1; j < N; j++) {
				// get previousXPrime
				previousXPrime[0] = xPrime[j-1];
				previousXPrime[1] = xDotPrime[j-1];
				previousXPrime[2] = xDDotPrime[j-1];
				previousXPrime[3] = xPrime[j];
				previousXPrime[4] = xDotPrime[j];
				previousXPrime[5] = xDDotPrime[j];
				previousXPrime[6] = xPrime[j+1];
				previousXPrime[7] = xDotPrime[j+1];
				previousXPrime[8] = xDDotPrime[j+1];
				
				B_xDot[j] = sloshingEquation.getValue(t, previousXPrime);
			}
			B_xDot[N] = 0;
			
			//compute C_x, C_xDot
			// compute xPrime, xDotPrime
			for (int k = 0; k < N; k++) {
				xPrime[k] = xPrevious[k] + deltaT / 2 * B_x[k];
				xDotPrime[k] = xDotPrevious[k] + deltaT / 2 * B_xDot[k];
			}
			// get C_x
			for (int k = 0; k < N; k++) {
				C_x[k] = xDotPrime[k];
			}
			// compute C_xDot
			C_xDot[0] = 0;
			for (int j = 1; j < N; j++) {
				// get previousXPrime
				previousXPrime[0] = xPrime[j-1];
				previousXPrime[1] = xDotPrime[j-1];
				previousXPrime[2] = xDDotPrime[j-1];
				previousXPrime[3] = xPrime[j];
				previousXPrime[4] = xDotPrime[j];
				previousXPrime[5] = xDDotPrime[j];
				previousXPrime[6] = xPrime[j+1];
				previousXPrime[7] = xDotPrime[j+1];
				previousXPrime[8] = xDDotPrime[j+1];
							
				C_xDot[j] = sloshingEquation.getValue(t, previousXPrime);
			}
			C_xDot[N] = 0;
			
			//compute D_x, D_xDot
			// compute xPrime, xDotPrime
			for (int k = 0; k < N; k++) {
				xPrime[k] = xPrevious[k] + deltaT * C_x[k];
				xDotPrime[k] = xDotPrevious[k] + deltaT * C_xDot[k];
			}
			// get D_x
			for (int k = 0; k < N; k++) {
				D_x[k] = xDotPrime[k];
			}
			// compute D_xDot
			D_xDot[0] = 0;
			for (int j = 1; j < N; j++) {
				// get previousXPrime
				previousXPrime[0] = xPrime[j-1];
				previousXPrime[1] = xDotPrime[j-1];
				previousXPrime[2] = xDDotPrime[j-1];
				previousXPrime[3] = xPrime[j];
				previousXPrime[4] = xDotPrime[j];
				previousXPrime[5] = xDDotPrime[j];
				previousXPrime[6] = xPrime[j+1];
				previousXPrime[7] = xDotPrime[j+1];
				previousXPrime[8] = xDDotPrime[j+1];
										
				D_xDot[j] = sloshingEquation.getValue(t, previousXPrime);
			}
			D_xDot[N] = 0;
			
			// compute current x, xDot
			for (int n = 0; n < N + 1; n++) {
				x[n] = xPrevious[n] + deltaT / 6 * (A_x[n] + 2*B_x[n] + 2*C_x[n] + D_x[n]);
				xDot[n] = xPrevious[n] + deltaT / 6 * (A_xDot[n] + 2*B_xDot[n] + 2*C_xDot[n] + D_xDot[n]);
			}
			
			//System.out.println("x:" + "  " + x[0] + "  " + x[1] + "  " + x[2] + "  " + x[3] + "  " + x[4]);
			//System.out.println("xDot:" + "  " + xDot[0] + "  " + xDot[1] + "  " + xDot[2] + "  " + xDot[3] + "  " + xDot[4]);
					
			
			// compute current xDDot
			xDDot[0] = 0;
			for (int j = 1; j < N; j++) {
				// get current x, xDot and previous xDDot
				previousX[0] = x[j-1];
				previousX[1] = xDot[j-1];
				previousX[2] = xDDotPrevious[j-1];
				previousX[3] = x[j];
				previousX[4] = x[j];
				previousX[5] = xDDotPrevious[j];
				previousX[6] = x[j+1];
				previousX[7] = x[j+1];
				previousX[8] = xDDotPrevious[j+1];
				
				xDDot[j] = sloshingEquation.getValue(t, previousX);
			}
			xDDot[N] = 0;
			
			// store all results into 2D array
			for (int j = 0; j < N + 1; j++) {
				allXResult[rowCounter][0] = allTs[rowCounter][0];
				allXResult[rowCounter][j * 3 + 1] = x[j];
				allXResult[rowCounter][j * 3 + 2] = xDot[j];
				allXResult[rowCounter][j * 3 + 3] = xDDot[j];
			}			
			
			etas[rowCounter][0] = allTs[rowCounter][0]; // value of current t
			// compute eta
			for (int j = 1; j < N + 1; j++) {
				double eta = sloshingEquation.getEta(x[j], x[j-1]);
				etas[rowCounter][j] = eta;
			}
			
			// output to csv files
			//toCSVfile(allXResult, "x.csv");
			//toCSVfile(etas, "etas.csv");
			
			// update
			rowCounter++;
			t += deltaT;
			for (int j = 0; j < N + 1; j++) {
				xPrevious[j] = x[j];
				xDotPrevious[j] = xDot[j];
				xDDotPrevious[j] = xDDot[j];
			}
			
		} // end of while
	} // end of main
	
	private static void toCSVfile(double[][] value, String filename) 
	{
		//Delimiter used in CSV file
				final String COMMA_DELIMITER = ","; // this was supposed t be declared as 'public static final'
				final String NEW_LINE_SEPARATOR = "\n"; // this was supposed t be declared as 'public static final'
				
			    //Code to save output in csv file
			    FileWriter fileWriter = null;
			    
			    //Code to save output in csv file
			    
				try
				{
					
					fileWriter = new FileWriter(filename);
					for (int i=0; i<value.length; i++) {
						for (int j=0; j<value[i].length; j++) {
//							System.out.printf(" %11f ", value[i][j]);
							fileWriter.append(String.valueOf(value[i][j]));
							fileWriter.append(COMMA_DELIMITER);
						}
						fileWriter.append(NEW_LINE_SEPARATOR);
//						System.out.println();
					}
				}
				catch (Exception e)
		    	{
		    		System.out.println("Error in CsvFileWriter!");
		    		e.printStackTrace();
		    	}
		    		
		    	finally 
		    	{
		    		try 
		    		{
		    			fileWriter.flush();
		    			fileWriter.close();
		    		} 
		    		catch (IOException e) 
		    		{
		    				System.out.println("Error while flushing/closing fileWriter!");
		                    e.printStackTrace();
		    		}
		    			
			    } 
				
			}
} // end of class
