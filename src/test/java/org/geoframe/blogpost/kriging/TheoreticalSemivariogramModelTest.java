package org.geoframe.blogpost.kriging;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;

import org.geoframe.blogpost.kriging.variogram.theoretical.model.*;

public class TheoreticalSemivariogramModelTest {
	double dist1 = 5.0;
	double dist2 = 25.0;

	double sill = 10.0;
	double range = 15.0;
	double nug = 1.0;
	private static final double EPSILON = 1e-2;

	// Helper method to check that the gradient array is valid.
	private void assertValidGradient(Model model) {
		double[] gradient = model.computeGradient();
		assertNotNull(gradient, "Il gradient non deve essere null");
		assertEquals(3, gradient.length, "Il gradient deve contenere 3 elementi");
	}

	@Test
	public void testLinearSemivariance() {

		Linear linear = new Linear(dist1, sill, range, nug);
		double computed = linear.computeSemivariance();
		assertEquals(4.33, computed, EPSILON, "Linear semivariance did not match expected value");
		linear = new Linear(dist2, sill, range, nug);
		computed = linear.computeSemivariance();
		assertEquals(11, computed, EPSILON, "Linear semivariance did not match expected value");
		assertValidGradient(linear);
	}

//    @Test
//    public void testBesselSemivariance() {
//        double dist = 5.0;
//        double sill = 10.0;
//        double range = 15.0;
//        double nug = 1.0;
//        Bessel bessel = new Bessel(dist, sill, range, nug);
//        double computed = bessel.computeSemivariance();
//        // Instead of assertNotEquals, use assertFalse with Double.compare
//        assertFalse("Bessel semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(bessel);
//    }

//    @Test
//    public void testCircularSemivariance() {
//        double dist = 5.0;
//        double sill = 10.0;
//        double range = 15.0;
//        double nug = 1.0;
//        Circular circular = new Circular(dist, sill, range, nug);
//        double computed = circular.computeSemivariance();
//        assertFalse("Circular semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(circular);
//    }

	@Test
	public void testExponentialSemivariance() {
		Exponential exponential = new Exponential(dist1, sill, range, nug);
		double computed = exponential.computeSemivariance();
		assertEquals(3.834, computed, EPSILON, "Exponential semivariance did not match expected value");
		exponential = new Exponential(dist2, sill, range, nug);
		computed = exponential.computeSemivariance();
		assertEquals(9.111, computed, EPSILON, "Exponential semivariance did not match expected value");
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Exponential semivariance should be computed");
		assertValidGradient(exponential);
	}

	@Test
	public void testGaussianSemivariance() {
		Gaussian gaussian = new Gaussian(dist1, sill, range, nug);
		double computed = gaussian.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Gaussian semivariance should be computed");
		assertEquals(2.051, computed, EPSILON, "Gausian semivariance did not match expected value");

		gaussian = new Gaussian(dist2, sill, range, nug);
		computed = gaussian.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Gaussian semivariance should be computed");
		assertEquals(10.378, computed, EPSILON, "Gaussian semivariance did not match expected value");

		assertValidGradient(gaussian);
	}

//    @Test
//    public void testHoleSemivariance() {
//        double dist = 5.0;
//        double sill = 8.0;
//        double range = 12.0;
//        double nug = 0.0;
//        Hole hole = new Hole(dist, sill, range, nug);
//        double computed = hole.computeSemivariance();
//        assertFalse("Hole semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(hole);
//    }

	@Test
	public void testLogarithmicSemivariance() {
		double dist = 5.0;
		double sill = 7.0;
		double range = 10.0;
		double nug = 0.5;
		Logarithmic logarithmic = new Logarithmic(dist, sill, range, nug);
		double computed = logarithmic.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Logarithmic semivariance should be computed");
		assertValidGradient(logarithmic);
	}

//    @Test
//    public void testPhentaspheriacalSemivariance() {
//        double dist = 5.0;
//        double sill = 9.0;
//        double range = 14.0;
//        double nug = 1.0;
//        Pentaspherical pentaspheriacal = new Pentaspherical(dist, sill, range, nug);
//        double computed = pentaspheriacal.computeSemivariance();
//        assertFalse("Phentaspheriacal semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(pentaspheriacal);
//    }

//    @Test
//    public void testPeriodicSemivariance() {
//        double dist = 5.0;
//        double sill = 11.0;
//        double range = 18.0;
//        double nug = 0.5;
//        Periodic periodic = new Periodic(dist, sill, range, nug);
//        double computed = periodic.computeSemivariance();
//        assertFalse("Periodic semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(periodic);
//    }

	@Test
	public void testPowerSemivariance() {
		double dist = 5.0;
		double sill = 10.0;
		double range = 15.0;
		double nug = 0.0;
		Power power = new Power(dist, sill, range, nug);
		double computed = power.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Power semivariance should be computed");
		assertValidGradient(power);
	}

	@Test
	public void testSphericalSemivariance() {
		Spherical spherical = new Spherical(dist1, sill, range, nug);
		double computed = spherical.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0);
		assertEquals(5.8148, computed, EPSILON, "Spherical semivariance should be computed");
		spherical = new Spherical(dist2, sill, range, nug);
		computed = spherical.computeSemivariance();
		assertFalse(Double.compare(Double.MAX_VALUE, computed) == 0, "Spherical semivariance should be computed");
		assertEquals(11.0, computed, EPSILON, "Spherical semivariance did not match expected value");

		assertValidGradient(spherical);
	}

//    @Test
//    public void testSplineSemivariance() {
//        double dist = 5.0;
//        double sill = 10.0;
//        double range = 15.0;
//        double nug = 0.5;
//        Spline spline = new Spline(dist, sill, range, nug);
//        double computed = spline.computeSemivariance();
//        assertFalse("Spline semivariance should be computed", Double.compare(Double.MAX_VALUE, computed) == 0);
//        assertValidGradient(spline);
//    }

}
