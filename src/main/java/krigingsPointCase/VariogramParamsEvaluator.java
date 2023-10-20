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

	@In
	public double[] x;

	@In
	public double[] y;

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

	@Out
	public double relError;

	@Out
	public String outSemivariogramType;

	@Execute
	public void proces() {
		try {

			ArrayList<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();

			if (expVar != null) {
				expVar.process();

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
			} else if (x != null && y != null) {
				for (int i = 0; i < x.length; i++) {
					double distance = x[i];
					double variance = y[i];
					WeightedObservedPoint point = new WeightedObservedPoint(1.0, distance, variance);
					points.add(point);
				}
			}

			String[] variogramType;
			if (pSemivariogramType != null) {
				variogramType = new String[] { pSemivariogramType };
			} else {
				variogramType = Utility.availableTheorethicalVariogra;
			}

			performEvaluation(variogramType, points);

			// System.out.println("calculated value " + nugget + " sill " + sill + " range "
			// + range);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			isFitGood = false;
		}
	}

	private void performEvaluation(String[] variogramType, ArrayList<WeightedObservedPoint> points) {
		relError = Double.MAX_VALUE;
		for (int j = 0; j < variogramType.length; j++) {
			try {

				VariogramFunction variogramFunction = new VariogramFunction(variogramType[j]);
				VariogramFitter fitter = new VariogramFitter(variogramFunction);

				ArrayList<WeightedObservedPoint> filteresPoints = variogramFunction.filterPoint(points);
				double coeffs[] = fitter.fit(filteresPoints);

				double distance = 0;
				int count = 0;
				for (int i = 0; i < points.size(); i++) {
					double actualValue = points.get(i).getY();
					if (actualValue != 0) {
						distance = distance + (SimpleModelFactory
								.createModel(variogramType[j], points.get(i).getX(), coeffs[0], coeffs[1], coeffs[2])
								.computeSemivariance() - actualValue) / actualValue;
						count = count + 1;
					}
				}
				double tmpError = Math.abs(distance / points.size());
				// System.out.println(variogramType[j]+" errore "+tmpError);
				if (j == 0 || tmpError < relError) {
					relError = Math.abs(distance / points.size());
					sill = coeffs[0];
					range = coeffs[1];
					nugget = coeffs[2];
					rmse = fitter.getRMS();
					isFitGood = relError < TRESHOLD;
					outSemivariogramType = variogramType[j];
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				System.out.println(e.getMessage());
			}
		}
	}
}
