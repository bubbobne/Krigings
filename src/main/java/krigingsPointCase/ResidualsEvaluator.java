package krigingsPointCase;

import flanagan.analysis.Regression;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

public class ResidualsEvaluator {
	@Description("Switch for detrended mode.")
	@In
	public boolean doDetrended;

	@Description("Switch for detrended mode.")
	@In
	public int regressionOrder = 1;

	@In
	public double[] zStations;
	@In
	public double[] hStations;

	@Out
	public double trend_intercept;

	@Out
	public double trend_coefficient;

	@In
	@Out
	public double[] hResiduals;

	@Out
	public boolean isPValueOk = false;

	@Execute
	public void process() {
		hResiduals = hStations;
		trend_intercept = 0;
		trend_coefficient = 0;
		isPValueOk = false;

		try {
			if (doDetrended && zStations.length > 2) {

				Regression r = new Regression();
				r = new Regression(zStations, hStations);
				r.polynomial(regressionOrder);
				/*
				 * If there is a trend for meteorological variables and elevation and it is
				 * statistically significant then the residuals from this linear trend are
				 * computed for each meteorological stations.
				 */
				if (r.getPvalues()[1] < 0.05) {
					isPValueOk = true;
					trend_intercept = r.getBestEstimates()[0];
					trend_coefficient = r.getBestEstimates()[1];
					hResiduals = r.getResiduals();
				} else {
					// it's set to true at each time step
					// set to 0 so the trend (line 330) is 0
					trend_intercept = 0;
					trend_coefficient = 0;
					hResiduals = hStations;
				}
			}
		} catch (Exception e) {
			// TODO: handle exception
			isPValueOk = false;
			hResiduals = hStations;
			trend_intercept = 0;
			trend_coefficient = 0;
		}

	}

}
