/* This file is part of JGrasstools (http://www.jgrasstools.org)
 * (C) HydroloGIS - www.hydrologis.com 
 * 
 * JGrasstools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package krigingsPointCase;

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.libs.exceptions.ModelsRuntimeException;
import org.hortonmachine.gears.libs.modules.HMConstants;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.math.matrixes.ColumnVector;
import org.hortonmachine.hmachine.i18n.HortonMessageHandler;
import org.hortonmachine.hmachine.modules.statistics.kriging.variogram.VariogramFunction;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;

import curvefitter.VariogramFitter;
import experimentalVariogram.ExperimentalVariogram;
import linearSistemSolver.SimpleLinearSystemSolverFactory;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
import theoreticalVariogram.TheoreticalVariogram;

@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello, Marialaura Bancheri & Francesco Serafin")
@Keywords("Kriging, Hydrology")
@Label("")
@Name("kriging")
@Status()
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class Krigings extends HMModel {

	@Description("The .shp of the measurement point, containing the position of the stations.")
	@In
	public SimpleFeatureCollection inStations = null;

	@Description("The field of the vector of stations, defining the id.")
	@In
	public String fStationsid = null;

	@Description("The field of the vector of stations, defining the elevation.")
	@In
	public String fStationsZ = null;

	@Description("The type of theoretical semivariogram: exponential, gaussian, spherical, pentaspherical"
			+ "linear, circular, bessel, periodic, hole, logaritmic, power, spline")
	@In
	public String pSemivariogramType = null;

	@Description("The HM with the measured data to be interpolated.")
	@In
	public HashMap<Integer, double[]> inData = null;

	@Description("The vector of the points in which the data have to be interpolated.")
	@In
	public SimpleFeatureCollection inInterpolate = null;

	@Description("The field of the interpolated vector points, defining the id.")
	@In
	public String fInterpolateid = null;

	@Description("The field of the interpolated vector points, defining the elevation.")
	@In
	public String fPointZ = null;

	@Description("The progress monitor.")
	@In
	public IHMProgressMonitor pm = new LogProgressMonitor();

	@Description("Include zeros in computations (default is true).")
	@In
	public boolean doIncludeZero = true;

	@Description("The range if the models runs with the gaussian variogram.")
	@In
	public double range;

	@Description("The sill if the models runs with the gaussian variogram.")
	@In
	public double sill;

	@Description("Is the nugget if the models runs with the gaussian variogram.")
	@In
	public double nugget;

	@Description("In the case of kriging with neighbor, maxdist is the maximum distance "
			+ "within the algorithm has to consider the stations")
	@In
	public double maxdist;

	@Description("In the case of kriging with neighbor, inNumCloserStations is the number "
			+ "of stations the algorithm has to consider")
	@In
	public int inNumCloserStations;

	@Description("Switch for detrended mode.")
	@In
	public boolean doDetrended;

	@Description("Degree of polynomial regression, default is 1")
	@In
	private final static int REGRESSION_ORDER = 1;

	@Description("Type of linear system solver")
	@In
	public String linearSystemSolverType = "default";

	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outData = null;

	@Description("Specified cutoff")
	@In
	public double cutoffInput;

	@Description("Number of bins to consider in the anlysis")
	@In
	public int cutoffDivide;
	@Description("Distances input path")
	@In
	public String inHValuesPath;

	/** transform to log. */
	@In
	public boolean doLogarithmic = false;

	/** transform to log. */
	@In
	public boolean boundedToZero = false;
	private static final double TOLL = 1.0d * 10E-8;

	private int step = 0;

	private HortonMessageHandler msg = HortonMessageHandler.getInstance();

	/** The id of the cosidered station */
	int id;

	@In
	public double pAGlobal;
	@In
	public double pSGlobal;

	@In
	public double pNugGlobal;

	@In
	public double pAGlobalDeTrended;

	@In
	public double pSGlobalDeTrended;

	@In
	public double pNugGlobalDeTrended;
	@In
	public String fileNoValue = "-9999";

	@In
	public String tStart = null;

	@In
	public int tTimeStep = 60;

	private double aGlobal;

	private double sGlobal;

	private double nugGlobal;

	private double aGlobalDeTrended;

	private double sGlobalDeTrended;

	private double nugGlobalDeTrended;

	private boolean noGlobalParams = false;

	/**
	 * Executing ordinary kriging.
	 * <p>
	 * <li>Verify if the parameters are correct.
	 * <li>Calculating the matrix of the covariance (a).
	 * <li>For each point to interpolated, evalutate the know term vector (b) and
	 * solve the system (a x)=b where x is the weight.
	 * </p>
	 *
	 * @throws Exception the exception
	 */

	@Execute
	public void executeKriging() throws Exception {

		verifyInput();
		nugget = 0;
		sill = 0;
		range = 0;
		LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;

		pointsToInterpolateId2Coordinates = getCoordinate(0, inInterpolate, fInterpolateid);

		Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();

		int j = 0;

		double[] result = new double[pointsToInterpolateId2Coordinates.size()];
		int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];

		/**
		 * StationsSelection is an external class that allows the selection of the
		 * stations involved in the study. It is possible to define if to include
		 * stations with zero values, station in a define neighborhood or within a max
		 * distance from the considered point.
		 */
		StationsSelection stations = new StationsSelection();
		stations.inStations = inStations;
		stations.doIncludezero = doIncludeZero;
		stations.maxdist = maxdist;
		stations.fStationsid = fStationsid;
		stations.fStationsZ = fStationsZ;
		stations.doLogarithmic = doLogarithmic;

		if (step == 0) {
			try {
				this.createDefaulParams(stations);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				pm.errorMessage(e.toString());
			}
		}
		step = step + 1;
		stations.inData = inData;
		stations.execute();

		double[] xStations = stations.xStationInitialSet;
		double[] yStations = stations.yStationInitialSet;
		double[] zStations = stations.zStationInitialSet;
		double[] hStations = stations.hStationInitialSet;

		int n1 = xStations.length - 1;
		double[] hResiduals = new double[n1];
		boolean areAllEquals = stations.areAllEquals;
		double trendCoeff = 0;
		double trendIntercept = 0;
		ResidualsEvaluator residualsEvaluator;
		if (!areAllEquals && n1 > 1) {
			residualsEvaluator = getResidualsEvaluator(Arrays.copyOfRange(zStations, 0, n1),
					Arrays.copyOfRange(hStations, 0, n1));
			hResiduals = residualsEvaluator.hResiduals;
			trendCoeff = residualsEvaluator.trend_coefficient;
			trendIntercept = residualsEvaluator.trend_intercept;
			if (!areAllEquals && n1 > 1) {
				setKrigingParams(stations, residualsEvaluator);
			}

		}
		while (idIterator.hasNext()) {
			double sum = 0.;
			id = idIterator.next();
			idArray[j] = id;
			Coordinate coordinate = (Coordinate) pointsToInterpolateId2Coordinates.get(id);
			if (maxdist > 0 || inNumCloserStations > 0) {
				if (inNumCloserStations > 0) {
					stations.inNumCloserStations = inNumCloserStations;
				}
				if (maxdist > 0) {
					stations.maxdist = maxdist;

				}
				stations.idx = coordinate.x;
				stations.idy = coordinate.y;
				stations.execute();
				xStations = stations.xStationInitialSet;
				yStations = stations.yStationInitialSet;
				zStations = stations.zStationInitialSet;
				hStations = stations.hStationInitialSet;
				n1 = xStations.length - 1;
			}

			xStations[n1] = coordinate.x;
			yStations[n1] = coordinate.y;
//			zStations[n1] = Double.NaN;

			if (n1 != 0) {

				if (!areAllEquals && n1 > 1) {

					pm.beginTask(msg.message("kriging.working"), pointsToInterpolateId2Coordinates.size());

					double h0 = 0.0;

					/*
					 * calculating the covariance matrix.
					 */
					double[][] covarianceMatrix = covMatrixCalculating(xStations, yStations, zStations, n1);

					double[] knownTerm = knownTermsCalculation(xStations, yStations, zStations, n1);

					/*
					 * solve the linear system, where the result is the weight
					 * (moltiplicativeFactor).
					 */
					ColumnVector solution = SimpleLinearSystemSolverFactory.solve(knownTerm, covarianceMatrix,
							linearSystemSolverType);

					double[] moltiplicativeFactor = solution.copyValues1D();

					for (int k = 0; k < n1; k++) {
						h0 = h0 + moltiplicativeFactor[k] * hResiduals[k];

						// sum is computed to check that
						// the sum of all the weights is 1
						sum = sum + moltiplicativeFactor[k];

					}

					// System.out.println(doDetrended);
					if (doDetrended) {
						double trend = coordinate.z * trendCoeff + trendIntercept;
						h0 = h0 + trend;
					}
					result[j] = h0;
					j++;

					if (Math.abs(sum - 1) >= TOLL) {
						throw new ModelsRuntimeException("Error in the coffeicients calculation",
								this.getClass().getSimpleName());
					}
					pm.worked(1);
				} else if (n1 == 1 || areAllEquals) {

					double tmp = hResiduals[0];
					pm.message(msg.message("kriging.setequalsvalue"));
					pm.beginTask(msg.message("kriging.working"), pointsToInterpolateId2Coordinates.size());
					result[j] = tmp;
					j++;
					n1 = 0;
					pm.worked(1);

				}

				pm.done();

			} else {

				pm.errorMessage("No value for this time step");

				double[] value = inData.values().iterator().next();
				result[j] = value[0];
				j++;

			}
			if (!isNovalue(result[j - 1])) {
				if (doLogarithmic) {
					result[j - 1] = Utility.getInverseLog(result[j - 1]);
				}
				if (boundedToZero && result[j - 1] < 0) {
					result[j - 1] = 0;
				}
			}
		}
		storeResult(result, idArray);

	}

	private ResidualsEvaluator getResidualsEvaluator(double[] zStations, double[] hStations) {
		ResidualsEvaluator residualsEvaluator = new ResidualsEvaluator();
		residualsEvaluator.doDetrended = this.doDetrended;
		residualsEvaluator.hStations = hStations;
		residualsEvaluator.zStations = zStations;
		residualsEvaluator.regressionOrder = this.REGRESSION_ORDER;
		residualsEvaluator.process();
		return residualsEvaluator;
	}

	private void setKrigingParams(StationsSelection stations, ResidualsEvaluator residualsEvaluator) throws Exception {

		int n1 = stations.hStationInitialSet.length - 1;
		double[] hResiduals = residualsEvaluator.hResiduals;
		if (this.inNumCloserStations != 0) {
			stations.inNumCloserStations = 0;
			stations.execute();
			n1 = stations.hStationInitialSet.length - 1;
			residualsEvaluator = getResidualsEvaluator(Arrays.copyOfRange(stations.zStationInitialSet, 0, n1),
					Arrays.copyOfRange(stations.hStationInitialSet, 0, n1));
			residualsEvaluator.process();

		}
		int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);
		hResiduals = residualsEvaluator.hResiduals;
		if (nugget == 0 && sill == 0 && range == 0) {
			VariogramParamsEvaluator variogramParamsEvaluator = new VariogramParamsEvaluator();
			variogramParamsEvaluator.expVar = getExperimentalVariogram(hResiduals, idStations);
			variogramParamsEvaluator.pSemivariogramType = this.pSemivariogramType;
			variogramParamsEvaluator.proces();

			if (residualsEvaluator.isPValueOk == false) {
				if ((nugget >= 0 && sill > 0 && range > 0 && variogramParamsEvaluator.isFitGood) || noGlobalParams) {
					nugget = variogramParamsEvaluator.nugget;
					sill = variogramParamsEvaluator.sill;
					range = variogramParamsEvaluator.range;
				} else {
					nugget = nugGlobal;
					sill = sGlobal;
					range = aGlobal;

				}
			} else if (residualsEvaluator.isPValueOk) {
				if ((nugget >= 0 && sill > 0 && range > 0 && variogramParamsEvaluator.isFitGood) || noGlobalParams) {
					nugget = variogramParamsEvaluator.nugget;
					sill = variogramParamsEvaluator.sill;
					range = variogramParamsEvaluator.range;
				} else {
					nugget = nugGlobalDeTrended;
					sill = sGlobalDeTrended;
					range = aGlobalDeTrended;

				}
			}
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
		if (inNumCloserStations > 0) {
			expVariogram.inNumCloserStations = this.inNumCloserStations;
		}
		if (cutoffDivide > 0) {
			expVariogram.Cutoff_divide = this.cutoffDivide;
		}
		if (cutoffInput > 0) {
			expVariogram.Cutoffinput = this.cutoffInput;
		}
		return expVariogram;
	}

	/**
	 * Verify the input of the model.
	 */
	private void verifyInput() {
		if (inData == null || inStations == null) {
			throw new NullPointerException(msg.message("kriging.stationProblem"));
		}
		if (doDetrended) {
			if (fPointZ == null || fStationsZ == null) {
				throw new NullPointerException("z field not found");
			}

			int ff = inStations.getSchema().indexOf(fStationsZ);
			int ff2 = inInterpolate.getSchema().indexOf(fPointZ);

			if (ff < 0 || ff2 < 0) {
				throw new NullPointerException("check if the z field name is correct");

			}

		}

		noGlobalParams = (inHValuesPath == null || inHValuesPath.isEmpty());
		if (nugget == 0 && range == 0 && nugget == 0 && noGlobalParams) {
			pm.message("Warning: nugget,range and sill will be evaluate only at each timestep");
		}

	}

	private void createDefaulParams(StationsSelection stations) throws Exception {
		nugGlobal = pNugGlobal;
		aGlobal = pAGlobal;
		sGlobal = pSGlobal;
		nugGlobalDeTrended = pNugGlobalDeTrended;
		sGlobalDeTrended = pSGlobalDeTrended;
		aGlobalDeTrended = pAGlobalDeTrended;

		if (!noGlobalParams && nugGlobal == 0 && aGlobal == 0 && sGlobal == 0) {
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

							if (this.inNumCloserStations != 0) {
								stations.inNumCloserStations = 0;

							}
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
					}

					nRows = nRows + 1;

				}

				readH.close();

				if (distance.length != variance.length) {
					throw new IllegalArgumentException(msg.message("kriging.defaultMode"));
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
				nugGlobal = coeffs[2];
				sGlobal = coeffs[0];
				aGlobal = coeffs[1];

				pm.message("Global value for nugget: " + nugGlobal + " sill:" + sGlobal + " range: " + aGlobal);

				if (doDetrended) {
					if (distanceDeTrended.length != varianceDeTrended.length) {
						throw new IllegalArgumentException(msg.message("kriging.defaultMode"));
					}
					if (nRowsDeTrended > 0) {

						points = new ArrayList<WeightedObservedPoint>();

						for (int i = 0; i < distanceDeTrended.length; i++) {
							double x = distanceDeTrended[i] / nRowsDeTrended;
							double y = varianceDeTrended[i] / nRowsDeTrended;
							WeightedObservedPoint point = new WeightedObservedPoint(1.0, x, y);
							points.add(point);
						}

						coeffs = fitter.fit(points);
						nugGlobalDeTrended = coeffs[2];
						sGlobalDeTrended = coeffs[0];
						aGlobalDeTrended = coeffs[1];

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

	/**
	 * Round.
	 *
	 * @param value  is the value of the variable considered
	 * @param places the places to consider after the comma
	 * @return the double value of the variable rounded
	 */
	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		long factor = (long) Math.pow(10, places);
		value = value * factor;
		long tmp = Math.round(value);
		return (double) tmp / factor;
	}

	/**
	 * Extract the coordinate of a FeatureCollection in a HashMap with an ID as a
	 * key.
	 *
	 * @param nStaz      the number of the stations
	 * @param collection is the collection of the considered points
	 * @param idField    the field containing the id of the stations
	 * @return the coordinate of the station
	 * @throws Exception if a field of elevation isn't the same of the collection
	 */
	private LinkedHashMap<Integer, Coordinate> getCoordinate(int nStaz, SimpleFeatureCollection collection,
			String idField) throws Exception {
		LinkedHashMap<Integer, Coordinate> id2CoordinatesMcovarianceMatrix = new LinkedHashMap<Integer, Coordinate>();
		FeatureIterator<SimpleFeature> iterator = collection.features();
		Coordinate coordinate = null;
		try {
			while (iterator.hasNext()) {
				SimpleFeature feature = iterator.next();
				int name = ((Number) feature.getAttribute(idField)).intValue();
				coordinate = ((Geometry) feature.getDefaultGeometry()).getCentroid().getCoordinate();
				double z = 0;
				if (fPointZ != null) {
					try {
						z = ((Number) feature.getAttribute(fPointZ)).doubleValue();
					} catch (NullPointerException e) {
						pm.errorMessage(msg.message("kriging.noPointZ"));
						throw new Exception(msg.message("kriging.noPointZ"));
					}
				}
				coordinate.z = z;
				id2CoordinatesMcovarianceMatrix.put(name, coordinate);
			}
		} finally {
			iterator.close();
		}

		return id2CoordinatesMcovarianceMatrix;
	}

	/**
	 * Covariance matrix calculation.
	 *
	 * @param x the x coordinates.
	 * @param y the y coordinates.
	 * @param z the z coordinates.
	 * @param n the number of the stations points.
	 * @return the double[][] matrix with the covariance
	 */
	private double[][] covMatrixCalculating(double[] x, double[] y, double[] z, int n) {

		double[][] covarianceMatrix = new double[n + 1][n + 1];

		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				double rx = x[i] - x[j];
				double ry = y[i] - y[j];
				double rz = z[i] - z[j];

				covarianceMatrix[j][i] = variogram(nugget, range, sill, rx, ry, rz);
				covarianceMatrix[i][j] = variogram(nugget, range, sill, rx, ry, rz);

			}
		}

		for (int i = 0; i < n; i++) {
			covarianceMatrix[i][n] = 1.0;
			covarianceMatrix[n][i] = 1.0;

		}
		covarianceMatrix[n][n] = 0;
		return covarianceMatrix;

	}

	/**
	 * Known terms calculation.
	 *
	 * @param x the x coordinates.
	 * @param y the y coordinates.
	 * @param z the z coordinates.
	 * @param n the number of the stations points.
	 * @return the double[] vector of the known terms
	 */
	private double[] knownTermsCalculation(double[] x, double[] y, double[] z, int n) {

		// known terms vector
		double[] gamma = new double[n + 1];

		for (int i = 0; i < n; i++) {
			double rx = x[i] - x[n];
			double ry = y[i] - y[n];
			double rz = z[i] - z[n];
			gamma[i] = variogram(nugget, range, sill, rx, ry, rz);
		}

		gamma[n] = 1.0;
		return gamma;

	}

	/**
	 * Variogram.
	 *
	 * @param nug   is the nugget
	 * @param range is the range
	 * @param sill  is the sill
	 * @param rx    is the x distance
	 * @param ry    is the y distance
	 * @param rz    is the z distance
	 * @return the double value of the variance
	 */
	private double variogram(double nug, double range, double sill, double rx, double ry, double rz) {
		if (isNovalue(rz)) {
			rz = 0;
		}
		double h2 = Math.sqrt(rx * rx + rz * rz + ry * ry);
		double vgmResult;

		if (h2 != 0) {
			TheoreticalVariogram vgm = new TheoreticalVariogram();
			vgmResult = vgm.calculateVGM(pSemivariogramType, h2, sill, range, nug);
		} else {
			vgmResult = 0;
		}
		return vgmResult;
	}

	/**
	 * Store the result in a HashMcovarianceMatrix (if the mode is 0 or 1).
	 *
	 * @param result the result
	 * @param id     the associated id of the calculating points.
	 * @throws SchemaException the schema exception
	 */
	private void storeResult(double[] result, int[] id) throws SchemaException {
		outData = new HashMap<Integer, double[]>();
		for (int i = 0; i < result.length; i++) {
			outData.put(id[i], new double[] { result[i] });
		}
	}

}
