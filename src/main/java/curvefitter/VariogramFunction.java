package curvefitter;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;

import theoreticalVariogram.Model;
import theoreticalVariogram.SimpleModelFactory;

public class VariogramFunction implements ParametricUnivariateFunction {

	String type;

	public VariogramFunction(String type) {
		this.type = type;
	}

	@Override
	public double value(double x, double... parameters) {
		// TODO: check if there is 3 parameters.

		Model variogram = SimpleModelFactory.createModel(type, x, parameters[0], parameters[1], parameters[2]);

		return variogram.computeSemivariance();
	}

	@Override
	public double[] gradient(double x, double... parameters) {
		
		Model variogram = SimpleModelFactory.createModel(type, x, parameters[0], parameters[1], parameters[2]);

		
		return variogram.computeGradient();

	}

}
