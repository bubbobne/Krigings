package org.geoframe.blogpost.kriging;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import org.geoframe.blogpost.kriging.linearsystemsolver.SimpleLinearSystemSolverFactory;
import org.geoframe.blogpost.kriging.pointcase.KrigingPointCase;
import org.geoframe.blogpost.kriging.variogram.theoretical.VariogramParameters;
import org.hortonmachine.gears.utils.math.matrixes.ColumnVector;

public class KrigingFunctionTest {

	@Test
	public void testCovarianceMatrixCalculationUsingReflection() throws Exception {
		// result from R
		// [,1] [,2] [,3]
		// [1,] 0.0000000 0.6321206 1
		// [2,] 0.6321206 0.0000000 1
		// [3,] 1.0000000 1.0000000 0

		// Synthetic data for two stations
		// Create variogram parameters with an exponential model:
		// Variogram model: gamma(h) = nugget + sill * (1 - exp(-h/range))
		// For h = 1, semivariance = 1 - exp(-1) ≈ 0.6321, with nugget = 0, sill = 1,
		// range = 1.

		double[] x = { 0.0, 1.0, 0.0 };
		// Create variogram parameters with an exponential model:
		// Variogram model: gamma(h) = nugget + sill * (1 - exp(-h/range))
		// For h = 1, semivariance = 1 - exp(-1) ≈ 0.6321, with nugget = 0, sill = 1,
		// range = 1.

		double[] y = { 0.0, 0.0, 1.0 };
		double[] z = { 0.0, 0.0, 0.0 };
		int n = 2; // number of stations

		// Instantiate your Kriging class and set the variogram parameters.
		Kriging kriging = new KrigingPointCase();
		VariogramParameters vp = new VariogramParameters.Builder("exponential", 0.0, 1.0, 1.0).setLocal(false)
				.setTrend(false).build();

		// Get the private field 'variogramParameters' from the Kriging class
		Field vpField = Kriging.class.getDeclaredField("variogramParameters");
		// Allow access to the private field
		vpField.setAccessible(true);
		// Set the field's value to our vp object
		vpField.set(kriging, vp);
		// Use reflection to access the non-public covMatrixCalculating method.
		Method method = Kriging.class.getDeclaredMethod("covMatrixCalculating", double[].class, double[].class,
				double[].class, int.class);
		method.setAccessible(true); // Allow access to the non-public method

		// Invoke the method
		double[][] covMatrix = (double[][]) method.invoke(kriging, x, y, z, n);

		// Expected semivariance for distance 1.
		double expectedSemivariance = 1.0 - Math.exp(-1.0);

		// The covariance matrix should be of size (n+1) x (n+1) = 3x3.
		// For i, j in {0, 1}:
		// - covMatrix[0][0] = 0 (distance 0)
		// - covMatrix[0][1] = expectedSemivariance (distance 1)
		// - covMatrix[1][0] = expectedSemivariance (distance 1)
		// - covMatrix[1][1] = 0 (distance 0)
		assertEquals(0.0, covMatrix[0][0], 1e-6);
		assertEquals(expectedSemivariance, covMatrix[0][1], 1e-6);
		assertEquals(expectedSemivariance, covMatrix[1][0], 1e-6);
		assertEquals(0.0, covMatrix[1][1], 1e-6);

		// The extra row/column (index 2) for the interpolation point:
		// - covMatrix[i][2] = 1.0 and covMatrix[2][i] = 1.0 for i = 0,1
		// - covMatrix[2][2] = 0
		assertEquals(1.0, covMatrix[0][2], 1e-6);
		assertEquals(1.0, covMatrix[1][2], 1e-6);
		assertEquals(1.0, covMatrix[2][0], 1e-6);
		assertEquals(1.0, covMatrix[2][1], 1e-6);
		assertEquals(0.0, covMatrix[2][2], 1e-6);
		Method method2 = Kriging.class.getDeclaredMethod("knownTermsCalculation", double[].class, double[].class,
				double[].class, int.class);
		method2.setAccessible(true);
		double[] knownTerm = (double[]) method2.invoke(kriging, x, y, z, n);
		assertEquals(0.6321, knownTerm[0], 1e-3);
		assertEquals(0.7568, knownTerm[1], 1e-3);
	}

	@Test
	public void testSolveSystem() throws Exception {
		double[] x = { 0.0, 0.0, 10.0, 10.0, 15.0 };
		double[] y = { 0.0, 10.0, 0.0, 10.0, 15.0 };
		double[] z = { 0.0, 0.0, 0.0, 0.0, 0.0 };
		int n = 4; // number of stations

		// Instantiate your Kriging class and set the variogram parameters.
		Kriging kriging = new KrigingPointCase();
		VariogramParameters vp = new VariogramParameters.Builder("exponential", 0.0, 1.0, 1.0).setLocal(false)
				.setTrend(false).build();

		// Get the private field 'variogramParameters' from the Kriging class
		Field vpField = Kriging.class.getDeclaredField("variogramParameters");
		// Allow access to the private field
		vpField.setAccessible(true);
		// Set the field's value to our vp object
		vpField.set(kriging, vp);
		// Use reflection to access the non-public covMatrixCalculating method.
		Method method = Kriging.class.getDeclaredMethod("covMatrixCalculating", double[].class, double[].class,
				double[].class, int.class);
		method.setAccessible(true); // Allow access to the non-public method

		// Invoke the method
		double[][] covMatrix = (double[][]) method.invoke(kriging, x, y, z, n);
		Method method2 = Kriging.class.getDeclaredMethod("knownTermsCalculation", double[].class, double[].class,
				double[].class, int.class);
		method2.setAccessible(true);
		double[] knownTerm = (double[]) method2.invoke(kriging, x, y, z, n);
		ColumnVector solution = SimpleLinearSystemSolverFactory.solve(knownTerm, covMatrix, "default");
		solution.print();
		assertEquals(0.25, solution.at(0), 1e-3);
		assertEquals(0.25, solution.at(1), 1e-3);
		assertEquals(0.25, solution.at(2), 1e-3);
		assertEquals(0.25, solution.at(3), 1e-3);
		assertEquals(0.24917, solution.at(4), 1e-3);
	}

}
