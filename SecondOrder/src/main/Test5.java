
public class Test5 {
	public static void main(String[] args) {

		FlyingBall ball = new FlyingBall();
		SecondDerivative acceleration = (SecondDerivative) ball;

		double t0 = 0;
		double tmax = 10;
		double y0 = 0;
		double h = 0.5;
		double v0 = 50;
		double a0 = acceleration.getValue(t0, y0, v0);
		IntegratorRKN_new integratorRKN = new IntegratorRKN_new(acceleration, t0, y0, v0, h, a0);

		System.out.println("TEST5, equation of motion");
		System.out.println(" t,s     y,m   v,m/s  a,m/s2");
		System.out.println("----------------------------");

		double t = t0, y = y0, v = v0, a = a0;
		double[][] RKN_result = new double[21][10];
		int indexRow = 0, indexCol = 0;

		RKN_result[0][0] = t;
		RKN_result[0][1] = y;
		RKN_result[0][2] = v;
		RKN_result[0][3] = a;

		do {
			// First row, initial states
			if (t == 0) {
				for (int i = 0; i < RKN_result[0].length; i++) {
					System.out.printf(" %f ", RKN_result[0][i]); // Default values are 0
				}
				System.out.println();
				integratorRKN.updateT0(); // update t to next time interval
				indexRow++; // go to next row
			}

			// Calculate all x, dx, d2x for each time interval
			// i is hard coded for now
			for (int i = 0; i < 3; i++) {
				integratorRKN.step(t, y, v, a);
				t = integratorRKN.t0;
				y = integratorRKN.y0;
				v = integratorRKN.dy0;
				a = integratorRKN.d2y0;

				// Store results in the array
				if (i == 0) {
					RKN_result[indexRow][indexCol] = t; // Store t once in one row
					System.out.printf(" %f ", t);
				}
				RKN_result[indexRow][++indexCol] = y;
				RKN_result[indexRow][++indexCol] = v;
				RKN_result[indexRow][++indexCol] = a;

				System.out.printf(" %f %f %f ", y, v, a);
			}

			// Update t0 and go to next row of the array, which is also next time interval
			integratorRKN.updateT0();
			++indexRow;
			indexCol = 0;

			System.out.println();

		} while (t < tmax);
		
		System.out.println();
		System.out.println();
		
		for (int i=0; i<RKN_result.length; i++) {
			for (int j=0; j<RKN_result[i].length; j++) {
				System.out.printf(" %f ", RKN_result[i][j]);
			}
			System.out.println();
		}
	}
}
