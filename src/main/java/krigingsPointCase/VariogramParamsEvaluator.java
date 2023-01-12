package krigingsPointCase;

import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.hortonmachine.gears.libs.modules.HMConstants;

import curvefitter.VariogramFitter;
import curvefitter.VariogramFunction;
import experimentalVariogram.ExperimentalVariogram;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import theoreticalVariogram.SimpleModelFactory;

public class VariogramParamsEvaluator {

	private static final double TRESHOLD = 0.5;

	@In
	public String pSemivariogramType;

	@In
	public ExperimentalVariogram expVar;

	@Out
	public double sill;

	@Out
	public double nugget;

	@Out
	public double range;

	@Out
	public double rmse;

	@Out
	public boolean isFitGood;

	@Execute
	public void proces() {
		try {
			expVar.process();

			VariogramFitter fitter = new VariogramFitter(new VariogramFunction(pSemivariogramType));
			ArrayList<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();
			Set<Entry<Integer, double[]>> entrySet = expVar.outDistances.entrySet();
			for (Entry<Integer, double[]> entry : entrySet) {
				Integer ID = entry.getKey();
				double x = expVar.outDistances.get(ID)[0];
				double y = expVar.outExperimentalVariogram.get(ID)[0];
				if (y != HMConstants.doubleNovalue) {
					WeightedObservedPoint point = new WeightedObservedPoint(1.0, x, y);
					points.add(point);
				}
			}

			double coeffs[] = fitter.fit(points);

			sill = coeffs[0];
			range = coeffs[1];
			nugget = coeffs[2];
			rmse = fitter.getRMS();
			double distance = 0;
			int count = 0;
			for (int i = 0; i < points.size(); i++) {
				double actualValue = points.get(i).getY();
				if (actualValue != 0) {
					distance = distance + (SimpleModelFactory
							.createModel(pSemivariogramType, points.get(i).getX(), sill, range, nugget)
							.computeSemivariance() - actualValue) / actualValue;
					count = count + 1;
				}
			}
			isFitGood = distance / points.size() < TRESHOLD;
			// System.out.println("calculated value " + nugget + " sill " + sill + " range "
			// + range);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			isFitGood = false;
		}
	}

}
