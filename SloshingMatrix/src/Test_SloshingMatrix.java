import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.linear.*;

public class Test_SloshingMatrix {
	public static Equation_SloshingMatrix sloshingEquation = new Equation_SloshingMatrix();
	public static int N = sloshingEquation.N; // use 4 for now  N = sloshingEquation.N;
	public static double[][] matrixData = new double[N-1][N-1]; // Matrix M

	public static void main(String[] args) {
		
		double t = 0;
		double tMax = sloshingEquation.T_MAX;
		double deltaT = sloshingEquation.DELTA_T;		
		double m = sloshingEquation.m;
		double x = 0;
		double xDot = 0;		
		double intermidiateXDDot = 0;     //value of F(x,xdot)
		double eta = 0;

		SecondDerivative_SloshingMatrix secondDerivative = (SecondDerivative_SloshingMatrix) sloshingEquation;
		IntegratorRKN_SloshingMatrix integratorRKN = new IntegratorRKN_SloshingMatrix (secondDerivative, t, x, xDot, deltaT, intermidiateXDDot);
		
		double[][] RKN_result = new double[(int)(tMax/deltaT)+1][3*(N+1)+1]; //[(int)(tMax/deltaT)+1][3*(N+1)+1]
		double[][] etas = new double[(int)(tMax/deltaT)+1][N+1]; // [(int)(tMax/deltaT)+1][N+1]to store ETA
		int indexRow = 0, indexCol = 3; // indexCol does not start from 0 since we need to leave columns for t and x0, xdot0, xdotdot0
		double [] previousX = new double[9]; // For new equation, we are not going to use all the nine values. 
			
		double[] fData = new double[N-1]; // to build vector F
		//double[][] matrixData = new double[N-1][N-1]; // Matrix M
		double[] xDDots = new double[N-1]; // to store xDDot for each time interval
	
		// t0, first row, initial states, all "0"
		for (int i = 0; i < RKN_result[0].length; i++) {
			RKN_result[0][i] = 0;
		}
		for (int i = 0; i < etas[0].length; i++) {
			etas[0][i] = 0;
		}

		integratorRKN.updateT0(); // update t to next time interval
		indexRow++; // go to the 2nd row	
		double xPrevious; // to store previous x value which will be used to calculate ETA

		for (int i = 1; i <= (int)(tMax / deltaT); i++) {			
			// Calculate all x, dx, d2x for each time interval
			// i is hard coded for now, should be i <= tMax / deltaT
			
			t = integratorRKN.t0;
			RKN_result[indexRow][0] = t;
			etas[indexRow][0] = t;
			
			// get x, xdot
			for (int j = 1; j <= (N-1); j++) { 
				
				// previous x value for calculating ETA later
				if (j == 1) {
					xPrevious = 0; // for the first ETA, previous x value is always 0
				} else {
					xPrevious = integratorRKN.y0; 
				}				
				
				previousX = getPreviousXValues(RKN_result, i, j);
				integratorRKN.step(t, x, xDot, intermidiateXDDot, previousX);
				t = integratorRKN.t0;
				x = integratorRKN.y0;
				xDot = integratorRKN.dy0;
				intermidiateXDDot = integratorRKN.d2y0;
				
				// Store x, xdot into final result 2D array array
				RKN_result[indexRow][++indexCol] = x;
				RKN_result[indexRow][++indexCol] = xDot;
				indexCol ++; // leave one column for xddot
				//RKN_result[indexRow][++indexCol] = a;
				
				// store F value
				fData[j-1] = intermidiateXDDot;
				
				// calculate and store ETA
				//eta = sloshingEquation.getEta(x, xPrevious);
			    eta = sloshingEquation.getEta(x, xPrevious);
				etas[indexRow][j] = eta;
				
				// get all 9 values for next inner loop
				previousX = getPreviousXValues(RKN_result, i, j);

				//solve new equations to get xDDot				
				// fill in matrix M
				fillInMatrix (previousX, m, j);
			}
						
			xDDots = solveLinearSystems (matrixData, fData);
			
			// store xDDot into final result 2D array
			for (int k = 0; k < xDDots.length; k++) {
				RKN_result[i][6 + k * 3] = xDDots[k];
			}
			
			// the nth (last) eta
			etas[i][N] = sloshingEquation.getEta(0, integratorRKN.dy0);
			
			// Update t0 and go to next row of the array, which is also next time interval
			integratorRKN.updateT0();
			++indexRow;
			indexCol = 3;

		} 
		
		/*
		System.out.println();
		System.out.println();
		
		System.out.println("TEST, equation of motion");
		System.out.println("     t,s           x0          x0dot       x0dotdot       x1          x1dot       x1dotdot       x2         x2dot        x2dotdot        x3        x3dot       x3dotdot         x4        x4dot        x4dotdot");
		System.out.println("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
		
		for (int i=0; i<RKN_result.length; i++) {
			for (int j=0; j<RKN_result[i].length; j++) {
				System.out.printf(" %11f ", RKN_result[i][j]);
			}
			System.out.println();
			//toCSVfile(RKN_result, "x.csv");
			//System.out.println("CSV file created successfully in project folder!");
		}
		*/
		toCSVfile(RKN_result, "x.csv");
		
		/*
		System.out.println("\n        t           eta1        eta2         eta3         eta4");
		System.out.println("------------------------------------------------------------");		
		for (int i=0; i<etas.length; i++) {
			for (int j=0; j<etas[i].length; j++) {
				System.out.printf(" %11f ", etas[i][j]);	
			}
			System.out.println();
			//toCSVfile(etas, "etas.csv");
			//System.out.println("CSV file created successfully in project folder!");
		}
		*/
		toCSVfile(etas, "etas.csv");
	}
	
	
	// Get the 9 values in previous state, x, xdot, xdotdot for i-1, i, i+1
	public static double[] getPreviousXValues (double[][] allValues, int currentRow, int currentNumOfX) {
		double[] xValues = new double[9];
		int previousRow = currentRow - 1;
		if (currentNumOfX == 1) {
			xValues[0] = 0;
			xValues[1] = 0;
			xValues[2] = 0;
		} else {
			xValues[0] = allValues[previousRow][1+(currentNumOfX-2)*3];
		    xValues[1] = allValues[previousRow][2+(currentNumOfX-2)*3];
		    xValues[2] = allValues[previousRow][3+(currentNumOfX-2)*3];
		}
		
		xValues[3] = allValues[previousRow][1+(currentNumOfX-1)*3];
		xValues[4] = allValues[previousRow][2+(currentNumOfX-1)*3];
		xValues[5] = allValues[previousRow][3+(currentNumOfX-1)*3];
		
		if (currentNumOfX == N) {
			xValues[6] = 0;
			xValues[7] = 0;
			xValues[8] = 0;
		} else {
			xValues[6] = allValues[previousRow][1+(currentNumOfX)*3];
		    xValues[7] = allValues[previousRow][2+(currentNumOfX)*3];
		    xValues[8] = allValues[previousRow][3+(currentNumOfX)*3];
		}
		
		return xValues;
	}
	
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
	// x (solution.getEntry(0)), y (solution.getEntry(1)), and z (solution.getEntry(2)) that solve the system. 
	public static double[] solveLinearSystems (double[][] matrixData, double[] constantData) {
		
		RealMatrix coefficients = new Array2DRowRealMatrix(matrixData, false); // matrix M
		DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
		RealVector constants = new ArrayRealVector(constantData, false); // vector F
		RealVector solution = solver.solve(constants); 
		
		// convert RealVector to array
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
	
}
