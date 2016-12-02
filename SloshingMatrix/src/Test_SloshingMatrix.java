import org.apache.commons.math3.linear.*;

public class Test_SloshingMatrix {
	public static int N = 4;
	
	public static void main(String[] args) {

		Equation_SloshingMatrix sloshingEquation = new Equation_SloshingMatrix();

		SecondDerivative_SloshingMatrix secondDerivative = (SecondDerivative_SloshingMatrix) sloshingEquation;
		double[][] RKN_result = new double[20][3*(N+1)+1]; //[][3*(N+1)+1], use 20 for testing
		int indexRow = 0, indexCol = 3; // indexCol does not start from 0 since we need to leave columns for t and x0, xdot0, xdotdot0
		double [] previousX = new double[9]; // For new equation, we are going to use less previous values. This will be modified later.
		
		double[] fData = new double[N-1]; // to build vector F
		double[][] matrixData = new double[N-1][N-1]; // Matrix M
		double[] xDDot = new double[N-1]; // to store xDDot for each time interval
		
		double t0 = 0;
		double tMax = sloshingEquation.T_MAX;
		double deltaT = sloshingEquation.DELTA_T;		
		double m = sloshingEquation.m;
		double y0 = 0;
		double v0 = 0;		//dot
		double a0 = 0;     //value of F
		IntegratorRKN_SloshingMatrix integratorRKN = new IntegratorRKN_SloshingMatrix (secondDerivative, t0, y0, v0, deltaT, a0);
		
		System.out.println("TEST, equation of motion");
		System.out.println(" t,s           x       xdot     xdotdot");
		System.out.println("-----------------------------------------");
		
		// t0, first row, initial states, all "0"
		for (int i = 0; i < RKN_result[0].length; i++) {
			RKN_result[0][i] = 0;
			//System.out.printf(" %f ", RKN_result[0][i]);
		}
		//System.out.println();
		//System.out.println();
		integratorRKN.updateT0(); // update t to next time interval
		indexRow++; // go to the 2nd row

		double t = integratorRKN.t0, y = y0, v = v0, a = a0;		

		for (int i = 1; i < 20; i++) {			
			// Calculate all x, dx, d2x for each time interval
			// i is hard coded for now, should be i <= tMax/deltaT
			
			RKN_result[indexRow][0] = t;
			//System.out.printf(" %f ", t);
			
			// get x, xdot
			for (int j = 1; j <= (N-1); j++) { 
				integratorRKN.step(t, y, v, a, previousX);
				t = integratorRKN.t0;
				y = integratorRKN.y0;
				v = integratorRKN.dy0;
				a = integratorRKN.d2y0;
				
				// Store x, xdot in the array
				RKN_result[indexRow][++indexCol] = y;
				RKN_result[indexRow][++indexCol] = v;
				indexCol ++; // leave one column for xddot
				//RKN_result[indexRow][++indexCol] = a;
				
				fData[j-1] = a;
				
				// get all 9 values for next inner loop
				previousX = getPreviousXValues(RKN_result, i, j);

				//System.out.printf(" %f %f %f ", y, v, a);
				// print 9 value array
				//System.out.println();
				//for (int k = 0; k < 9; k++) {
				//	System.out.printf(" %f ", previousX[k]);
				//}
				//System.out.println();
				
				//solve new equations to get xDDot
				// form matrix M
				double gammaCurrent = sloshingEquation.getGammaCurrent (previousX); // Gamma i
				double gammaNext = sloshingEquation.getGammaNext (previousX); // Gamma i+1
				if (j == 1) {
					matrixData[j-1][0] = m + gammaCurrent + gammaNext;
					matrixData[j-1][1] = - gammaNext;
				} else if (j == (N-1)) {
					matrixData[j-1][j-2] = - gammaCurrent;
					matrixData[j-1][j-1] = m + gammaCurrent + gammaNext;
				} else {
					matrixData[j-1][j-2] = - gammaCurrent;
					matrixData[j-1][j-1] = m + gammaCurrent + gammaNext;
					matrixData[j-1][j] = - gammaNext;
				}				

			}
						
			xDDot = solveLinearSystems (matrixData, fData);
			// store xDDot into result 2D array
			for (int k = 0; k < xDDot.length; k++) {
				RKN_result[i][6 + k * 3] = xDDot[k];
			}

			
			// Update t0 and go to next row of the array, which is also next time interval
			integratorRKN.updateT0();
			++indexRow;
			indexCol = 3;

			//System.out.println();

		} 
		
		//System.out.println();
		//System.out.println();
		
		//for (int i=0; i<RKN_result.length; i++) {
		//	for (int j=0; j<RKN_result[i].length; j++) {
		//		System.out.printf(" %f ", RKN_result[i][j]);
		//	}
		//	System.out.println();
		//}

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
	
}
