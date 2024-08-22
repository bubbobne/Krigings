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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.hortonmachine.gears.libs.exceptions.ModelsRuntimeException;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.math.matrixes.ColumnVector;
import org.hortonmachine.hmachine.i18n.HortonMessageHandler;
import org.locationtech.jts.geom.Coordinate;
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
import theoreticalVariogram.VariogramParameters;

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
	@Description("The hashmap withe the interpolated results")
	@Out
	public HashMap<Integer, double[]> outData = null;

	@Description("The hashmap withe the parameter for each time step: 0=>nugget, 1=>sill, 2=> range, 3 => is local (0 if used local parameter 1 if use global, 4 => is detrended (0 no trend, 1 with trend)")
	@Out
	public HashMap<Integer, double[]> outVariogramParams = null;
	@Description("The Experimental Distances.")
	@Out
	public HashMap<Integer, double[]> outDistances;

	@Description("The Experimental Variogram.")
	@Out
	public HashMap<Integer, double[]> outExperimentalVariogram;
	@Description("The numbers of pairs at certain lag.")
	@Out
	public HashMap<Integer, double[]> outNumberPairsPerBin;

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

	private boolean noGlobalParams = false;
	@In
	public String globalDetrendedVariogramType;
	@In
	public String globalVariogramType;

	private boolean noEvaluationGlobalParams;

	private VariogramParameters variogramParameters;

	private VariogramParameters vpGlobal;

	private VariogramParameters vpGlobalTrend;

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

		if (step == 0) {
			try {
				verifyInput();
				this.createDefaulParams();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				pm.errorMessage(e.toString());
			}
		}
//		nugget = 0;
//		sill = 0;
//		range = 0;
		LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;

		pointsToInterpolateId2Coordinates = Utility.getCoordinate(inInterpolate, fInterpolateid, fPointZ, pm, msg);

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
		stations.maxdist = maxdist;
		stations.fStationsid = fStationsid;
		stations.fStationsZ = fStationsZ;
		stations.doLogarithmic = doLogarithmic;
		stations.doIncludezero = doIncludeZero;
		stations.inNumCloserStations = 0;
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
			variogramParameters = setKrigingParams(stations, residualsEvaluator);

		} else {
			variogramParameters = new VariogramParameters();
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
				stations.doIncludezero = doIncludeZero;
				stations.doLogarithmic = doLogarithmic;
				stations.inData = inData;

				stations.execute();
				xStations = stations.xStationInitialSet;
				yStations = stations.yStationInitialSet;
				zStations = stations.zStationInitialSet;
				hStations = stations.hStationInitialSet;
				n1 = xStations.length - 1;
				residualsEvaluator = getResidualsEvaluator(Arrays.copyOfRange(zStations, 0, n1),
						Arrays.copyOfRange(hStations, 0, n1));
				hResiduals = residualsEvaluator.hResiduals;
				trendCoeff = residualsEvaluator.trend_coefficient;
				trendIntercept = residualsEvaluator.trend_intercept;
			}

			xStations[n1] = coordinate.x;
			yStations[n1] = coordinate.y;
			// zStations[n1] = Double.NaN;

			if (n1 != 0) {

				if (!areAllEquals && n1 > 1) {

					// pm.beginTask(msg.message("kriging.working"),
					// pointsToInterpolateId2Coordinates.size());

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
					// pm.worked(1);
				} else if (n1 == 1 || areAllEquals) {

					double tmp = hResiduals[0];
					// pm.message(msg.message("kriging.setequalsvalue"));
					// pm.beginTask(msg.message("kriging.working"),
					// pointsToInterpolateId2Coordinates.size());
					result[j] = tmp;
					j++;
					n1 = 0;
					// pm.worked(1);

				}

				// pm.done();

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
		residualsEvaluator.regressionOrder = Krigings.REGRESSION_ORDER;
		residualsEvaluator.process();
		return residualsEvaluator;
	}

	private VariogramParameters setKrigingParams(StationsSelection stations, ResidualsEvaluator residualsEvaluator)
			throws Exception {

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
		double isLocal = 0.0;
		double isTrend = 0.0;
		int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);
		hResiduals = residualsEvaluator.hResiduals;
		VariogramParameters vp = null;
		if (nugget == 0 && sill == 0 && range == 0) {
			VariogramParamsEvaluator variogramParamsEvaluator = new VariogramParamsEvaluator();
			variogramParamsEvaluator.expVar = getExperimentalVariogram(hResiduals, idStations);
			variogramParamsEvaluator.pSemivariogramType = this.pSemivariogramType;
			variogramParamsEvaluator.proces();
			this.outDistances = variogramParamsEvaluator.expVar.outDistances;
			this.outNumberPairsPerBin = variogramParamsEvaluator.expVar.outNumberPairsPerBin;
			this.outExperimentalVariogram = variogramParamsEvaluator.expVar.outExperimentalVariogram;

			if (residualsEvaluator.isPValueOk == false) {
				if ((variogramParamsEvaluator.nugget >= 0 && variogramParamsEvaluator.sill > 0
						&& variogramParamsEvaluator.range > 0 && variogramParamsEvaluator.isFitGood)
						|| noGlobalParams) {
					vp = new VariogramParameters(variogramParamsEvaluator.outSemivariogramType,
							variogramParamsEvaluator.nugget, variogramParamsEvaluator.range,
							variogramParamsEvaluator.sill);
					vp.setIsLocal(true);

				} else if (doDetrended) {
					vp = vpGlobalTrend;
				} else {
					vp = vpGlobal;
				}
			} else if (residualsEvaluator.isPValueOk) {
				if ((variogramParamsEvaluator.nugget >= 0 && variogramParamsEvaluator.sill > 0
						&& variogramParamsEvaluator.range > 0 && variogramParamsEvaluator.isFitGood)
						|| noGlobalParams) {
					vp = new VariogramParameters(variogramParamsEvaluator.outSemivariogramType,
							variogramParamsEvaluator.nugget, variogramParamsEvaluator.range,
							variogramParamsEvaluator.sill);
					vp.setIsLocal(true);

				} else if (doDetrended) {
					vp = vpGlobalTrend;
				} else {
					vp = vpGlobal;
				}

			}
			vp.setIsTrend(doDetrended);
		} else {
			vp = new VariogramParameters(pSemivariogramType, nugget, range, sill);
			vp.setIsLocal(false);
			vp.setIsTrend(false);
		}

		return vp;

	}

	private ExperimentalVariogram getExperimentalVariogram(double[] hresiduals, int[] idArray) {
		ExperimentalVariogram expVariogram = Utility.getExperimentalVariogram(fStationsid, inStations, doIncludeZero,
				cutoffDivide, cutoffInput, 0);
		HashMap<Integer, double[]> tmpInData = new HashMap<Integer, double[]>();
		for (int i = 0; i < idArray.length; i++) {
			tmpInData.put(idArray[i], new double[] { hresiduals[i] });
		}
		expVariogram.inData = tmpInData;
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
		if ((nugget != 0 || sill != 0 || range != 0) && (pSemivariogramType == null || pSemivariogramType.isEmpty())) {
			throw new NullPointerException(
					"You provide incomplete fixed parameter check: nugget, sill, range or pSemivariogramType");
		}
	}

	private void createDefaulParams() throws Exception {

		noEvaluationGlobalParams = (inHValuesPath == null || inHValuesPath.isEmpty());

		boolean gloablZero = (pNugGlobal == 0 && pAGlobal == 0 && pSGlobal == 0);
		boolean gloablZeroDeTrended = (pNugGlobalDeTrended == 0 && pAGlobalDeTrended == 0 && pSGlobalDeTrended == 0);

		if (noEvaluationGlobalParams && !gloablZero && globalVariogramType == null) {
			throw new NullPointerException(
					"You provide incomplete fixed parameter for global: nugGlobal, sGlobal, rGlobal or globalVariogramType");
		}

		if (noEvaluationGlobalParams && !gloablZeroDeTrended && globalDetrendedVariogramType == null) {
			throw new NullPointerException(
					"You provide incomplete fixed parameter for global de trended: nugGlobal, sGlobal, rGlobal or globalVariogramType");
		}

		if (!doDetrended) {
			noGlobalParams = noEvaluationGlobalParams && gloablZero || globalVariogramType == null;
			noEvaluationGlobalParams = noEvaluationGlobalParams || (!gloablZero && globalVariogramType != null);
		} else {
			noGlobalParams = noEvaluationGlobalParams && (gloablZero || globalVariogramType == null
					|| gloablZeroDeTrended || globalDetrendedVariogramType == null);
			noEvaluationGlobalParams = noEvaluationGlobalParams || (!gloablZero && globalVariogramType != null
					&& !gloablZeroDeTrended && globalDetrendedVariogramType != null);
		}

		if (nugget == 0 && range == 0 && nugget == 0 && noGlobalParams) {
			pm.message("Warning: nugget,range and sill will be evaluate only at each timestep");
		}
		if (!noEvaluationGlobalParams) {

			GlobalParameterEvaluator gParam = new GlobalParameterEvaluator();
			gParam.inStations = inStations;
			gParam.doIncludeZero = doIncludeZero;
			gParam.maxdist = maxdist;
			gParam.fStationsid = fStationsid;
			gParam.fStationsZ = fStationsZ;
			gParam.doLogarithmic = doLogarithmic;
			gParam.cutoffDivide = cutoffDivide;
			gParam.inHValuesPath = this.inHValuesPath;
			gParam.fileNoValue = fileNoValue;
			gParam.tStart = tStart;
			gParam.tTimeStep = tTimeStep;
			gParam.pSemivariogramType = pSemivariogramType;
			gParam.doDetrended = doDetrended;
			gParam.cutoffInput = cutoffInput;
			gParam.execute();
			vpGlobal = new VariogramParameters(gParam.outGlobalVariogramType, gParam.pNugGlogal, gParam.pAGlobal,
					gParam.pSGlogal);
			vpGlobal.setIsLocal(false);
			vpGlobal.setIsTrend(false);
			vpGlobalTrend = new VariogramParameters(gParam.outGlobalDetrendedVariogramType, gParam.pNugGlogalDeTrended,
					gParam.pAGlobalDeTrended, gParam.pSGlogalDeTrended);
			vpGlobal.setIsLocal(false);
			vpGlobal.setIsTrend(false);
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
				covarianceMatrix[j][i] = TheoreticalVariogram.calculateVGMxyz(variogramParameters, x[i], y[i], z[i],
						x[j], y[j], z[j]);
				covarianceMatrix[i][j] = covarianceMatrix[j][i];
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
			gamma[i] = TheoreticalVariogram.calculateVGMxyz(variogramParameters, x[i], y[i], z[i], x[n], y[n], z[n]);
		}
		gamma[n] = 1.0;
		return gamma;

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
		outExperimentalVariogram = variogramParameters.toHashMap();
	}

}
