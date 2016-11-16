
public class Test_2DArray {
	public static int N = 3;
	
	public static void main(String[] args) {

		FlyingBall_2DArray ball = new FlyingBall_2DArray();
		SecondDerivative_2DArray acceleration = (SecondDerivative_2DArray) ball;
		double[][] RKN_result = new double[21][10];
		int indexRow = 0, indexCol = 0;
		double [] xValues = new double[9];
		
		double t0 = 0;
		double y0 = 0;
		double h = 0.5;
		double v0 = 50;		
		double a0 = 0;
		IntegratorRKN_2DArray integratorRKN = new IntegratorRKN_2DArray(acceleration, t0, y0, v0, h, a0);

		System.out.println("TEST, equation of motion");
		System.out.println(" t,s     y,m   v,m/s  a,m/s2");
		System.out.println("----------------------------");
		
		// First row, initial states, all "0"
		for (int i = 0; i < RKN_result[0].length; i++) {
			RKN_result[0][i] = 0;
			System.out.printf(" %f ", RKN_result[0][i]);
		}
		System.out.println();
		integratorRKN.updateT0(); // update t to next time interval
		indexRow++; // go to the 2nd row

		double t = t0, y = y0, v = v0, a = a0;		

		for (int i = 1; i < 20; i++) {			
			// Calculate all x, dx, d2x for each time interval
			// i is hard coded for now
			for (int j = 1; j <= 3; j++) {
				integratorRKN.step(t, y, v, a, xValues);
				t = integratorRKN.t0;
				y = integratorRKN.y0;
				v = integratorRKN.dy0;
				a = integratorRKN.d2y0;

				// Store results in the array
				if (j == 1) {
					RKN_result[indexRow][indexCol] = t; // Store t once in one row
					System.out.printf(" %f ", t);
				}
				RKN_result[indexRow][++indexCol] = y;
				RKN_result[indexRow][++indexCol] = v;
				RKN_result[indexRow][++indexCol] = a;
				
				// get all 9 values for next inner loop
				xValues = getPreviousXValues(RKN_result, i, j);

				System.out.printf(" %f %f %f ", y, v, a);
				
				// print 9 value array to check if getting correct values
				System.out.println();
				for (int k = 0; k < 9; k++) {
					System.out.printf(" %f ", xValues[k]);
				}
				System.out.println();
			}

			// Update t0 and go to next row of the array, which is also next time interval
			integratorRKN.updateT0();
			++indexRow;
			indexCol = 0;

			System.out.println();

		} 
		
		System.out.println();
		System.out.println();
		
		// Print final array
		for (int i=0; i<RKN_result.length; i++) {
			for (int j=0; j<RKN_result[i].length; j++) {
				System.out.printf(" %f ", RKN_result[i][j]);
			}
			System.out.println();
		}
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
	
}
