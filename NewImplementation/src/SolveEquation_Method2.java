import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class SolveEquation_Method2 {
	static int N = 100;
	static EquationSloshingMatrix sloshingEquation = new EquationSloshingMatrix();
	public static double[][] matrixData = new double[N-1][N-1]; // Matrix M
	
	public static void main (String[] args) {
		double t = 0;
		double deltaT = 1.0 / 256;
		int T_MAX = 8;
		double m = sloshingEquation.m;
	
		double[] xPrevious = new double[N+1];
		double[] xDotPrevious = new double[N+1];
		double[] xDDotPrevious = new double[N+1];
		
		double[] x = new double[N+1];
		double[] xDot = new double[N+1];
		double[] xDDot = new double[N+1];
	
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
		
		double[] fData = new double[N-1]; // to build vector F
		
		while (t <= T_MAX) {
		//while (t <= 0.1) {
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
				
				fData[j-1] = sloshingEquation.getValue(t, previousX);				
				fillInMatrix(previousX, m, j);
			}			
			A_xDot[N] = 0;
			double[] temp = solveLinearSystems(matrixData, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				A_xDot[i+1] = temp[i];
			}
			
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
				
				fData[j-1] = sloshingEquation.getValue(t, previousXPrime);
				fillInMatrix(previousXPrime, m, j);
			}
			B_xDot[N] = 0;
			temp = solveLinearSystems(matrixData, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				B_xDot[i+1] = temp[i];
			}
			
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
							
				fData[j-1] = sloshingEquation.getValue(t, previousXPrime);
				fillInMatrix(previousXPrime, m, j);
			}
			C_xDot[N] = 0;
			temp = solveLinearSystems(matrixData, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				C_xDot[i+1] = temp[i];
			}
			
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
										
				fData[j-1] = sloshingEquation.getValue(t, previousXPrime);
				fillInMatrix(previousXPrime, m, j);
			}
			D_xDot[N] = 0;
			temp = solveLinearSystems(matrixData, fData);
			for (int i = 0; i < temp.length; i++) { // temp.length = N-1
				D_xDot[i+1] = temp[i];
			}
			
			// compute current x, xDot
			for (int n = 0; n < N + 1; n++) {
				x[n] = xPrevious[n] + deltaT / 6 * (A_x[n] + 2*B_x[n] + 2*C_x[n] + D_x[n]);
				xDot[n] = xPrevious[n] + deltaT / 6 * (A_xDot[n] + 2*B_xDot[n] + 2*C_xDot[n] + D_xDot[n]);
			}
			
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
				
				xDDot[j] = sloshingEquation.getValue(t, previousX); // previous DDot were not actually used
			}
			xDDot[N] = 0;
			
			// store all results into 2D array
			for (int j = 0; j < N + 1; j++) {
				allXResult[rowCounter][0] = allTs[rowCounter][0];
				allXResult[rowCounter][j * 3 + 1] = x[j];
				allXResult[rowCounter][j * 3 + 2] = xDot[j];
				allXResult[rowCounter][j * 3 + 3] = xDDot[j];
			}			
			
			// compute eta
			for (int j = 1; j < N + 1; j++) {
				double eta = sloshingEquation.getEta(x[j], x[j-1]);
				etas[rowCounter][0] = allTs[rowCounter][0];
				etas[rowCounter][j] = eta;
			}
			
			// output to csv files
			toCSVfile(allXResult, "x.csv");
			toCSVfile(etas, "etas.csv");
			
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
	
	public static void fillInMatrix (double[] previousX, double m, int index) {
		double gammaCurrent = sloshingEquation.getGammaCurrent (previousX); // Gamma i
		double gammaNext = sloshingEquation.getGammaNext (previousX); // Gamma i+1
		if (index == 1) {
			matrixData[index-1][0] = m + gammaCurrent + gammaNext;
			matrixData[index-1][1] = - gammaNext;
		} else if (index == (N-1)) {
			matrixData[index-1][index-2] = - gammaCurrent;
			matrixData[index-1][index-1] = m + gammaCurrent + gammaNext;
		} else {
			matrixData[index-1][index-2] = - gammaCurrent;
			matrixData[index-1][index-1] = m + gammaCurrent + gammaNext;
			matrixData[index-1][index] = - gammaNext;
		}
	}
	
	// Solve linear systems using matrix
	// The solution vector will contain values for 
	// x (solution.getEntry(0)), y (solution.getEntry(1)), and z (solution.getEntry(2))... that solve the system. 
	public static double[] solveLinearSystems (double[][] matrixData, double[] constantData) {
		
		RealMatrix coefficients = new Array2DRowRealMatrix(matrixData, false); // matrix M
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
