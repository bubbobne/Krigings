package KrigingsTests;
import static org.junit.Assert.assertEquals;

import java.lang.reflect.Field;
import java.lang.reflect.Method;

import org.junit.Test;
import org.geoframe.blogpost.kriging.pointcase.Kriging;
import org.geoframe.blogpost.kriging.variogram.theoretical.VariogramParameters;

public class TestKrigingFunction {

    	@Test
        public void testCovarianceMatrixCalculationUsingReflection() throws Exception {
            // Synthetic data for two stations
            // Create variogram parameters with an exponential model:
            // Variogram model: gamma(h) = nugget + sill * (1 - exp(-h/range))
            // For h = 1, semivariance = 1 - exp(-1) ≈ 0.6321, with nugget = 0, sill = 1, range = 1.

            double[] x = { 0.0, 1.0 };
            // Create variogram parameters with an exponential model:
            // Variogram model: gamma(h) = nugget + sill * (1 - exp(-h/range))
            // For h = 1, semivariance = 1 - exp(-1) ≈ 0.6321, with nugget = 0, sill = 1, range = 1.

            double[] y = { 0.0, 0.0 };
            double[] z = { 0.0, 0.0 };
            int n = 2; // number of stations

            // Instantiate your Kriging class and set the variogram parameters.
            Kriging kriging = new Kriging();
            VariogramParameters vp = new VariogramParameters.Builder("exponential", 0.0, 1.0, 1.0)
                    .setLocal(false)
                    .setTrend(false)
                    .build();

            // Get the private field 'variogramParameters' from the Kriging class
            Field vpField = Kriging.class.getDeclaredField("variogramParameters");
            // Allow access to the private field
            vpField.setAccessible(true);
            // Set the field's value to our vp object
            vpField.set(kriging, vp);
            // Use reflection to access the non-public covMatrixCalculating method.
            Method method = Kriging.class.getDeclaredMethod("covMatrixCalculating", double[].class, double[].class, double[].class, int.class);
            method.setAccessible(true);  // Allow access to the non-public method

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
        }
    
    
}
