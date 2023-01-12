package krigingsPointCase;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;

import curvefitter.VariogramFitter;
import curvefitter.VariogramFunction;
import experimentalVariogram.ExperimentalVariogram;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

public class GlobalParameterEvaluator {
	@Description("The .shp of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;
	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;

	@Description("The progress monitor.")
	@In
	public IHMProgressMonitor pm = new LogProgressMonitor();

	@Description("Include zeros in computations (default is true).")
	@In
	public boolean doIncludeZero = true;

	@Description("The type of theoretical semivariogram: exponential, gaussian, spherical, pentaspherical"
			+ "linear, circular, bessel, periodic, hole, logaritmic, power, spline")
	@In
	public String pSemivariogramType = null;
	@Description("Switch for detrended mode.")
	@In
	public boolean doDetrended;

	@Description("Specified cutoff")
	@In
	public double cutoffInput;

	@Description("Number of bins to consider in the anlysis")
	@In
	public int cutoffDivide;
	@Description("Distances input path")
	@In
	public String inHValuesPath;
	@In
	public String tStart = null;

	@In
	public int tTimeStep = 60;
	/** transform to log. */
	@In
	public boolean doLogarithmic = false;

	@Out
	public double pAGlobal;
	@Out

	public double pSGlogal;

	@Out
	public double pNugGlogal;

	@Out
	public double pAGlobalDeTrended;

	@Out
	public double pSGlogalDeTrended;

	@Out
	public double pNugGlogalDeTrended;
	@Description("In the case of kriging with neighbor, maxdist is the maximum distance "
			+ "within the algorithm has to consider the stations")
	@In
	public double maxdist;

	@In
	public String fileNoValue = "-9999";

	@Execute
	public void execute() {

		/**
		 * StationsSelection is an external class that allows the selection of the
		 * stations involved in the study. It is possible to define if to include
		 * stations with zero values, station in a define neighborhood or within a max
		 * distance from the considered point.
		 */

		if (doDetrended) {
			int ff = inStations.getSchema().indexOf(fStationsZ);

			if (ff < 0) {
				throw new NullPointerException("check if the z field name is correct");

			}

		}

		StationsSelection stations = new StationsSelection();
		stations.inStations = inStations;
		stations.doIncludezero = doIncludeZero;
		stations.maxdist = maxdist;
		stations.fStationsid = fStationsid;
		stations.fStationsZ = fStationsZ;
		stations.doLogarithmic = doLogarithmic;

		try {
			this.createDefaulParams(stations);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void createDefaulParams(StationsSelection stations) throws Exception {

		System.out.println("evaluate global params");
		double[] variance = new double[cutoffDivide];
		double[] distance = new double[cutoffDivide];

		double[] varianceDeTrended = new double[cutoffDivide];
		double[] distanceDeTrended = new double[cutoffDivide];
		OmsTimeSeriesIteratorReader readH = new OmsTimeSeriesIteratorReader();
		readH.file = this.inHValuesPath;
		readH.idfield = "ID";
		readH.fileNovalue = fileNoValue;
		readH.tStart = tStart;
		readH.tTimestep = tTimeStep;
		readH.initProcess();
		ExperimentalVariogram exp = getExperimentalVariogram();
		int nRows = 0;
		int nRowsDeTrended = 0;
		try {
			while (readH.doProcess) {
				readH.nextRecord();
				HashMap<Integer, double[]> h = readH.outData;

				exp.inData = h;
				exp.process();
				if (!exp.areAllEquals && exp.differents > 2) {
					HashMap<Integer, double[]> d = exp.outDistances;
					HashMap<Integer, double[]> v = exp.outExperimentalVariogram;
					int j = 0;
					for (Map.Entry<Integer, double[]> tt : d.entrySet()) {
						distance[j] = distance[j] + tt.getValue()[0];
						variance[j] = variance[j] + v.get(tt.getKey())[0];
						j = j + 1;
					}
					if (doDetrended) {

						stations.inData = h;
						stations.execute();
						int n1 = stations.hStationInitialSet.length - 1;
						int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);

						ResidualsEvaluator rEvaluator = new ResidualsEvaluator();
						rEvaluator.doDetrended = true;
						rEvaluator.hStations = Arrays.copyOfRange(stations.hStationInitialSet, 0, n1);
						rEvaluator.zStations = Arrays.copyOfRange(stations.zStationInitialSet, 0, n1);

						rEvaluator.process();
						if (rEvaluator.isPValueOk) {
							double[] hVal = rEvaluator.hResiduals;
							nRowsDeTrended = nRowsDeTrended + 1;
							HashMap<Integer, double[]> hRes = new HashMap<>();
							for (int i = 0; i < n1; i++) {
								hRes.put(stations.idStationInitialSet[i], new double[] { hVal[i] });
							}

							exp.inData = hRes;
							exp.process();
							if (exp.differents > 2) {
								d = exp.outDistances;
								v = exp.outExperimentalVariogram;
								j = 0;
								for (Map.Entry<Integer, double[]> tt : d.entrySet()) {
									distanceDeTrended[j] = distanceDeTrended[j] + tt.getValue()[0];
									varianceDeTrended[j] = varianceDeTrended[j] + v.get(tt.getKey())[0];
									j = j + 1;
								}
							}
						}
					}

					nRows = nRows + 1;

				}

			}

			readH.close();

			if (distance.length != variance.length) {
				throw new IllegalArgumentException();
			}

			VariogramFitter fitter = new VariogramFitter(new VariogramFunction(pSemivariogramType));
			ArrayList<WeightedObservedPoint> points = new ArrayList<WeightedObservedPoint>();

			for (int i = 0; i < distance.length; i++) {
				double x = distance[i] / nRows;
				double y = variance[i] / nRows;
				WeightedObservedPoint point = new WeightedObservedPoint(1.0, x, y);
				points.add(point);
			}

			double coeffs[] = fitter.fit(points);
			pNugGlogal = coeffs[2];
			pSGlogal = coeffs[0];
			pAGlobal = coeffs[1];

			pm.message("Global value for nugget: " + pNugGlogal + " sill:" + pSGlogal + " range: " + pAGlobal);

			if (doDetrended) {
				if (distanceDeTrended.length != varianceDeTrended.length) {
					throw new IllegalArgumentException();
				}

				points = new ArrayList<WeightedObservedPoint>();
				if (nRowsDeTrended > 0) {
					for (int i = 0; i < distanceDeTrended.length; i++) {
						double x = distanceDeTrended[i] / nRowsDeTrended;
						double y = varianceDeTrended[i] / nRowsDeTrended;
						WeightedObservedPoint point = new WeightedObservedPoint(1.0, x, y);
						points.add(point);
					}

					coeffs = fitter.fit(points);
					pNugGlogalDeTrended = coeffs[2];
					pSGlogalDeTrended = coeffs[0];
					pAGlobalDeTrended = coeffs[1];
					pm.message("Global value with TREND for nugget: " + pNugGlogalDeTrended + " sill:"
							+ pSGlogalDeTrended + " range: " + pAGlobalDeTrended);

				} else {
					pm.message("no trend has been found, so no parameters has been evauate");
				}
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private ExperimentalVariogram getExperimentalVariogram(double[] hresiduals, int[] idArray) {
		ExperimentalVariogram expVariogram = getExperimentalVariogram();
		HashMap<Integer, double[]> tmpInData = new HashMap<Integer, double[]>();
		for (int i = 0; i < idArray.length; i++) {
			tmpInData.put(idArray[i], new double[] { hresiduals[i] });
		}
		expVariogram.inData = tmpInData;
		return expVariogram;
	}

	private ExperimentalVariogram getExperimentalVariogram() {
		ExperimentalVariogram expVariogram = new ExperimentalVariogram();
		expVariogram.fStationsid = this.fStationsid;
		expVariogram.inStations = this.inStations;
		expVariogram.doIncludezero = this.doIncludeZero;

		expVariogram.inNumCloserStations = 0;

		if (cutoffDivide > 0) {
			expVariogram.Cutoff_divide = this.cutoffDivide;
		}
		if (cutoffInput > 0) {
			expVariogram.Cutoffinput = this.cutoffInput;
		}
		return expVariogram;
	}

}
