import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class SolveEquation_Method2 {
	static EquationSloshingMatrix sloshingEquation = new EquationSloshingMatrix();
	static int N = sloshingEquation.N;
	public static double[][] matrixM = new double[N-1][N-1]; // Matrix M
	
	public static void main (String[] args) {
		double t = 0;
		double deltaT = sloshingEquation.DELTA_T;
		int T_MAX = sloshingEquation.T_MAX;
		double m = sloshingEquation.m; // mass
	
		// for previous t
		double[] xPrevious = new double[N+1];
		double[] xDotPrevious = new double[N+1];
		double[] xDDotPrevious = new double[N+1];
		
		// for current t
		double[] x = new double[N+1];
		double[] xDot = new double[N+1];
		double[] xDDot = new double[N+1];
	
		// for RK4
		double[] xPrime = new double[N+1];
		double[] xDotPrime = new double[N+1];
		double[] xDDotPrime = new double[N+1];
	
		// coefficients
		double[] A_x = new double[N+1];
		double[] A_xDot = new double[N+1];
		double[] B_x = new double[N+1];
		double[] B_xDot = new double[N+1];
		double[] C_x = new double[N+1];
		double[] C_xDot = new double[N+1];
		double[] D_x = new double[N+1];
		double[] D_xDot = new double[N+1];
		
		//previousX[]= {x_i-1,   x_i-1_dot,     x_i-1_doubleDot,  
		//                     x_i,      xi_dot,          xi_doubleDot,   
		//                     x_i+1,   x_i+1_dot,    x_i+1_doubleDot}
		double[] previousX = new double[9];
		double[] previousXPrime = new double[9];
		
		double[][] allXResult = new double[(int)(T_MAX/deltaT)+1][3*(N+1)+1];
		double[][] etas = new double[(int)(T_MAX/deltaT)+1][N+1];
		double[][] allTs = new double[(int)(T_MAX/deltaT)+1][1];
		int rowCounter = 0;
		
		// vector F
		double[] fData = new double[N-1];
		
		while (t <= T_MAX) {
			allTs[rowCounter][0] = t;
			System.out.println("\nt: " + t);
			
			//******************************
			//* Calculate coefficient A *
			//******************************
			//get A_x
			// A_x[0], A_x[N] are always 0
			for (int i = 0; i < N + 1; i++) {
				A_x[i] = xDotPrevious[i];
			}			
			
			//compute A_xDot
			//A_xDot[0], A_xDot[N] are always 0
			for (int i = 1; i < N; i++) {
				// get previousX
				previousX = getPreviousX(i, xPrevious, xDotPrevious, xDDotPrevious);
				
				fData[i-1] = sloshingEquation.getValue(t, previousX);				
				fillInMatrixM(previousX, m, i);
			}			
			double[] temp = solveLinearSystems(matrixM, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				A_xDot[i+1] = temp[i];
			}

			//******************************
			//* Calculate coefficient B *
			//******************************
			// compute B_x, B_xDot
			// compute xPrime, xDotPrime, and get B_x
			for (int i = 0; i < N; i++) {
				xPrime[i] = xPrevious[i] + deltaT / 2 * A_x[i];
				xDotPrime[i] = xDotPrevious[i] + deltaT / 2 * A_xDot[i];
				B_x[i] = xDotPrime[i];
			}

			// compute B_xDot
			B_xDot[0] = 0;
			for (int i = 1; i < N; i++) {
				// get previousXPrime
				previousXPrime = getPreviousX(i, xPrime, xDotPrime, xDDotPrime);				
				
				fData[i-1] = sloshingEquation.getValue(t + deltaT / 2, previousXPrime);
				fillInMatrixM(previousXPrime, m, i);
			}			
			temp = solveLinearSystems(matrixM, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				B_xDot[i+1] = temp[i];
			}

			//******************************
			//* Calculate coefficient C *
			//******************************
			//compute C_x, C_xDot
			// compute xPrime, xDotPrime, and get C_x
			for (int i = 0; i < N; i++) {
				xPrime[i] = xPrevious[i] + deltaT / 2 * B_x[i];
				xDotPrime[i] = xDotPrevious[i] + deltaT / 2 * B_xDot[i];
				C_x[i] = xDotPrime[i];
			}
			
			//******************************
			//* Calculate coefficient C *
			//******************************
			// compute C_xDot
			C_xDot[0] = 0;
			for (int i = 1; i < N; i++) {
				// get previousXPrime
				previousXPrime = getPreviousX(i, xPrime, xDotPrime, xDDotPrime);
				
				fData[i-1] = sloshingEquation.getValue(t + deltaT / 2, previousXPrime);
				fillInMatrixM(previousXPrime, m, i);
			}
			temp = solveLinearSystems(matrixM, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				C_xDot[i+1] = temp[i];
			}

			//******************************
			//* Calculate coefficient D *
			//******************************
			//compute D_x, D_xDot
			// compute xPrime, xDotPrime, get D_x
			for (int i = 0; i < N; i++) {
				xPrime[i] = xPrevious[i] + deltaT * C_x[i];
				xDotPrime[i] = xDotPrevious[i] + deltaT * C_xDot[i];
				D_x[i] = xDotPrime[i];
			}

			// compute D_xDot
			D_xDot[0] = 0;
			for (int i = 1; i < N; i++) {
				// get previousXPrime
				previousXPrime = getPreviousX (i, xPrime, xDotPrime, xDDotPrime);
										
				fData[i-1] = sloshingEquation.getValue(t + deltaT, previousXPrime);
				fillInMatrixM(previousXPrime, m, i);
			}
			temp = solveLinearSystems(matrixM, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				D_xDot[i+1] = temp[i];
			}
			
			// compute current x, xDot
			for (int n = 0; n < N + 1; n++) {
				x[n] = xPrevious[n] + deltaT / 6 * (A_x[n] + 2*B_x[n] + 2*C_x[n] + D_x[n]);
				xDot[n] = xDotPrevious[n] + deltaT / 6 * (A_xDot[n] + 2*B_xDot[n] + 2*C_xDot[n] + D_xDot[n]);
			}

			//********************************
			//* Calculate current xDDot *
			//********************************
			xDDot[0] = 0;
			for (int i = 1; i < N; i++) {
				// get current x, xDot and previous xDDot
				previousX = getPreviousX(i, x, xDot, xDDotPrevious);
				
				xDDot[i] = sloshingEquation.getValue(t, previousX); // previous DDot were not actually used
			}
			xDDot[N] = 0;
			
			//********************************
			//* Calculate eta                  *
			//********************************
			for (int i = 1; i < N + 1; i++) {
				double eta = sloshingEquation.getEta(x[i], x[i-1]);
				etas[rowCounter][0] = allTs[rowCounter][0];
				etas[rowCounter][i] = eta;
				if (i==3)
				System.out.println("eta: " + eta);
			}
			
			//*************************************
			//* Save and output to csv files  *
			//*************************************
			// store all results into 2D array
			for (int i = 0; i < N + 1; i++) {
				allXResult[rowCounter][0] = allTs[rowCounter][0];
				allXResult[rowCounter][i * 3 + 1] = x[i];
				allXResult[rowCounter][i * 3 + 2] = xDot[i];
				allXResult[rowCounter][i * 3 + 3] = xDDot[i];
			}	
			// output to csv files
			toCSVfile(allXResult, "allX_newMethod.csv");
			toCSVfile(etas, "etas_newMethod.csv");
			
			// update
			rowCounter++;
			t += deltaT;
			for (int i = 0; i < N + 1; i++) {
				xPrevious[i] = x[i];
				xDotPrevious[i] = xDot[i];
				xDDotPrevious[i] = xDDot[i];
			}			
		} // end of while 		
	} // end of main
	
	//Get previous x, xdot, xddot
	//Return previousX[]= {x_i-1,   x_i-1_dot,     x_i-1_doubleDot,  
	//                                 x_i,      xi_dot,          xi_doubleDot,   
	//                                 x_i+1,   x_i+1_dot,    x_i+1_doubleDot}
	public static double[] getPreviousX (int i, double[] xPrevious, 
            													double[] xDotPrevious, double[] xDDotPrevious) {
		double[] previousX = new double[9];
		previousX[0] = xPrevious[i-1];
		previousX[1] = xDotPrevious[i-1];
		previousX[2] = xDDotPrevious[i-1];
		previousX[3] = xPrevious[i];
		previousX[4] = xDotPrevious[i];
		previousX[5] = xDDotPrevious[i];
		previousX[6] = xPrevious[i+1];
		previousX[7] = xDotPrevious[i+1];
		previousX[8] = xDDotPrevious[i+1];

		return previousX;
	}
	
	// Get matrix M for each i
	public static void fillInMatrixM (double[] previousX, double m, int index) {
		double gammaCurrent = sloshingEquation.getGammaCurrent (previousX); // Gamma i
		double gammaNext = sloshingEquation.getGammaNext (previousX); // Gamma i+1
		if (index == 1) {
			matrixM[index-1][0] = m + gammaCurrent + gammaNext;
			matrixM[index-1][1] = - gammaNext;
		} else if (index == (N-1)) {
			matrixM[index-1][index-2] = - gammaCurrent;
			matrixM[index-1][index-1] = m + gammaCurrent + gammaNext;
		} else {
			matrixM[index-1][index-2] = - gammaCurrent;
			matrixM[index-1][index-1] = m + gammaCurrent + gammaNext;
			matrixM[index-1][index] = - gammaNext;
		}
	}
	
	// Solve linear systems using matrix
	// The solution vector will contain values for 
	// x (solution.getEntry(0)), y (solution.getEntry(1)), and z (solution.getEntry(2))... that solve the system. 
	public static double[] solveLinearSystems (double[][] matrixM, double[] constantData) {
		
		RealMatrix coefficients = new Array2DRowRealMatrix(matrixM, false); // matrix M
		DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
		RealVector constants = new ArrayRealVector(constantData, false); // vector F
		RealVector solution = solver.solve(constants); 
		
		// convert RealVector solution to array
		int length = solution.getDimension();
		double[] xDDot = new double[length];
		
		for (int i = 0; i < length; i++) {
			xDDot [i] = solution.getEntry(i);
		}
		
		return xDDot;
	}
	
	// Save 2D array to csv file
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
							//System.out.printf(" %11f ", value[i][j]);
							fileWriter.append(String.valueOf(value[i][j]));
							fileWriter.append(COMMA_DELIMITER);
						}
						fileWriter.append(NEW_LINE_SEPARATOR);
						//System.out.println();
					}
					
					/* only output eta3
					for (int i=0; i<value.length; i++) {
						//for (int j=0; j<value[i].length; j++) {
							//System.out.printf(" %11f ", value[i][j]);
							fileWriter.append(String.valueOf(value[i][0]));
							fileWriter.append(COMMA_DELIMITER);
							fileWriter.append(String.valueOf(value[i][3]));
							fileWriter.append(COMMA_DELIMITER);
						//}
						fileWriter.append(NEW_LINE_SEPARATOR);
						//System.out.println();
					}
					*/
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
