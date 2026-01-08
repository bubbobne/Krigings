package org.geoframe.blogpost.kriging.variogram.theoretical.curvefitter;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.RealVector;

public class KrigingParamValidator implements ParameterValidator {
	final double[] lowerBounds_;

	public KrigingParamValidator(double[] lowerBounds) {
		if (lowerBounds.length != 3) {
			throw new DimensionMismatchException(lowerBounds.length, 3);
		}

		lowerBounds_ = lowerBounds;
	}

	@Override
	public RealVector validate(RealVector params) {
		if (params.getDimension() != 3) {
			throw new DimensionMismatchException(params.getDimension(), 3);
		}
		RealVector n = params.copy();
		for (int i = 0; i < 3; i++) {
			if (n.getEntry(i) < lowerBounds_[i]) {
				n.setEntry(i, lowerBounds_[i]);
			}
		}
		return n;
	}

}