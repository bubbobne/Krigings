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

	@Out
	public StationsSelection stations;

	@Description("In the case of kriging with neighbor, maxdist is the maximum distance "
			+ "within the algorithm has to consider the stations")
	@In
	public double maxdist;

	@In
	public String fileNoValue = "-9999";
	public String outGlobalVariogramType;
	public String outGlobalDetrendedVariogramType;

	@Execute
	public void execute() {
		verifyInput();

		stations = new StationsSelection();
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

	private void verifyInput() {
		if (doDetrended) {
			int ff = inStations.getSchema().indexOf(fStationsZ);

			if (ff < 0) {
				throw new NullPointerException("check if the z field name is correct");
			}
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
		ExperimentalVariogram exp = Utility.getExperimentalVariogram(fStationsid, inStations, doIncludeZero,
				cutoffDivide, cutoffInput, 0);
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

			for (int i = 0; i < distance.length; i++) {
				distance[i] = distance[i] / nRows;
				variance[i] = variance[i] / nRows;

			}
			VariogramParamsEvaluator vEvaluator = new VariogramParamsEvaluator();
			vEvaluator.pSemivariogramType = pSemivariogramType;
			vEvaluator.x = distance;
			vEvaluator.y = variance;
			vEvaluator.proces();

			pNugGlogal = vEvaluator.nugget;
			pSGlogal = vEvaluator.sill;
			pAGlobal = vEvaluator.range;
			outGlobalVariogramType = vEvaluator.outSemivariogramType;

			pm.message("Global value for nugget: " + pNugGlogal + " sill:" + pSGlogal + " range: " + pAGlobal
					+ "  semivariogram type:" + vEvaluator.outSemivariogramType);

			if (doDetrended) {
				if (distanceDeTrended.length != varianceDeTrended.length) {
					throw new IllegalArgumentException();
				}

				if (nRowsDeTrended > 0) {
					for (int i = 0; i < distanceDeTrended.length; i++) {
						distanceDeTrended[i] = distanceDeTrended[i] / nRowsDeTrended;
						varianceDeTrended[i] = varianceDeTrended[i] / nRowsDeTrended;
					}
					vEvaluator = new VariogramParamsEvaluator();
					vEvaluator.pSemivariogramType = pSemivariogramType;
					vEvaluator.x = distanceDeTrended;
					vEvaluator.y = varianceDeTrended;
					vEvaluator.proces();
					pNugGlogalDeTrended = vEvaluator.nugget;
					pSGlogalDeTrended = vEvaluator.sill;
					pAGlobalDeTrended = vEvaluator.range;
					outGlobalDetrendedVariogramType = vEvaluator.outSemivariogramType;

					pm.message("Global value with TREND for nugget: " + pNugGlogalDeTrended + " sill:"
							+ pSGlogalDeTrended + " range: " + pAGlobalDeTrended + "  semivariogram type:"
							+ vEvaluator.outSemivariogramType);

				} else {
					pm.message("no trend has been found, so no parameters has been evauate");
				}
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
