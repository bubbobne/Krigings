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
package org.geoframe.blogpost.kriging.pointcase;

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

import org.geoframe.blogpost.kriging.linearsystemsolver.SimpleLinearSystemSolverFactory;
import org.geoframe.blogpost.kriging.primarylocation.ResidualsEvaluator;
import org.geoframe.blogpost.kriging.primarylocation.StationProcessor;
import org.geoframe.blogpost.kriging.primarylocation.StationsSelection;
import org.geoframe.blogpost.kriging.utilities.Utility;
import org.geoframe.blogpost.kriging.variogram.theoretical.GlobalParameterEvaluator;
import org.geoframe.blogpost.kriging.variogram.theoretical.TheoreticalVariogram;
import org.geoframe.blogpost.kriging.variogram.theoretical.VariogramParameters;
import org.geoframe.blogpost.kriging.variogram.theoretical.VariogramParametersCalculator;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.hortonmachine.gears.libs.exceptions.ModelsRuntimeException;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.math.matrixes.ColumnVector;
import org.hortonmachine.gears.utils.math.matrixes.MatrixException;
import org.hortonmachine.hmachine.i18n.HortonMessageHandler;
import org.locationtech.jts.geom.Coordinate;

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

@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello, Marialaura Bancheri & Francesco Serafin")
@Keywords("Kriging, Hydrology")
@Label("")
@Name("kriging")
@Status()
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class Kriging extends HMModel {

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
	public final static int REGRESSION_ORDER = 1;

	@Description("Type of linear system solver")
	@In
	public String linearSystemSolverType = "math3";

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
	@Description("The hashmap with the interpolated results")
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
	public HashMap<Integer, double[]> outTheoreticalVariogram;
	@Description("The Experimental Variogram.")
	@In
	public HashMap<Integer, double[]> inTheoreticalVariogram;

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
	@In
	public String globalDetrendedVariogramType;
	@In
	public String globalVariogramType;

	private VariogramParameters variogramParameters;

	private VariogramParameters vpGlobal;

	private VariogramParameters vpGlobalTrend;

	/**
	 * Executes the Ordinary Kriging algorithm following these steps: 1. Input
	 * verification and variogram parameter setup. 2. Station selection based on
	 * distance and other criteria. 3. Calculation of the covariance matrix and the
	 * known terms vector. 4. Solving the linear system to obtain the weights. 5.
	 * Computing the interpolated value, applying the trend if detrending is
	 * enabled. 6. Final validation and storage of the results.
	 *
	 * @throws Exception if any error occurs during the kriging execution.
	 */

	@Execute
	public void executeKriging() throws Exception {
		// Initialization and verification
		VariogramParameters vp = initializeKrigingParameters();
		// Extract interpolation points from input feature collection
		LinkedHashMap<Integer, Coordinate> pointsToInterpolateId2Coordinates = null;
		pointsToInterpolateId2Coordinates = Utility.getCoordinate(inInterpolate, fInterpolateid, fPointZ, pm, msg);
		Set<Integer> pointsToInterpolateIdSet = pointsToInterpolateId2Coordinates.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSet.iterator();
		double[] result = new double[pointsToInterpolateId2Coordinates.size()];
		int[] idArray = new int[pointsToInterpolateId2Coordinates.size()];

		/**
		 * StationsSelection is an external class that allows the selection of the
		 * stations involved in the study. It is possible to define if to include
		 * stations with zero values, station in a define neighborhood or within a max
		 * distance from the considered point.
		 */
		StationsSelection stations = this.createAndConfigureStationSelection();
		determineVariogram(vp, stations);
		StationProcessor sp = new StationProcessor(stations);
		if (maxdist == 0 && inNumCloserStations == 0) {
			sp.updateForCoordinate(null, inData, 0, 0, doDetrended);
		}
		int j = 0;
		while (idIterator.hasNext()) {
			id = idIterator.next();
			idArray[j] = id;
			double[] xStations = null, yStations = null, zStations = null, hResiduals = null;
			Coordinate coordinate = (Coordinate) pointsToInterpolateId2Coordinates.get(id);
			if (maxdist > 0 || inNumCloserStations > 0) {
				sp.updateForCoordinate(coordinate, inData, inNumCloserStations, maxdist, doDetrended);
			}
			double trendIntercept = sp.getTrendIntercept();
			double trendCoeff = sp.getTrendCoeff();
			if (trendCoeff == 0 && trendIntercept == 0 && variogramParameters.getIsTrend()) {
				variogramParameters = vpGlobal;
			}

			int n1 = sp.getCount();
			boolean areAllEquals = sp.areAllEquals();
			double interpolatedValue = 0;

			// zStations[n1] = Double.NaN;
			if (n1 != 0) {

				if (!areAllEquals && n1 > 1) {

					interpolatedValue = intepolateValue(sp,coordinate);
					// pm.worked(1);
				} else if (n1 == 1 || areAllEquals) {

					interpolatedValue = hResiduals[0];
					// pm.message(msg.message("kriging.setequalsvalue"));
					// pm.beginTask(msg.message("kriging.working"),
					// pointsToInterpolateId2Coordinates.size());
					n1 = 0;
					// pm.worked(1);
				}

				// pm.done();

			} else {

				pm.errorMessage("No value for this time step");

				interpolatedValue = inData.values().iterator().next()[0];

			}
			result[j] = postProcessResult(interpolatedValue);
			j++;
		}
		storeResult(result, idArray);

	}

	private double intepolateValue(StationProcessor sp, Coordinate coordinate) throws MatrixException {
		double interpolatedValue;
		double[] xStations = sp.getXStations();
		double[] yStations = sp.getYStations();
		double[] zStations = sp.getZStations();
		double[] hResiduals = sp.getHResiduals();
		double trendCoeff = sp.getTrendCoeff();
		double trendIntercept = sp.getTrendIntercept();
		int n1 = sp.getCount();
		xStations[n1] = coordinate.x;
		yStations[n1] = coordinate.y;
		zStations[n1] = coordinate.z;
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
		double sum = 0;
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
		interpolatedValue = h0;

		if (Math.abs(sum - 1) >= TOLL) {
			throw new ModelsRuntimeException("Error in the coffeicients calculation", this.getClass().getSimpleName());
		}
		return interpolatedValue;
	}

	private double postProcessResult(double interpolatedValue) {
		if (!isNovalue(interpolatedValue)) {
			if (doLogarithmic) {
				interpolatedValue = Utility.getInverseLog(interpolatedValue);
			}
			if (boundedToZero && interpolatedValue < 0) {
				interpolatedValue = 0;
			}
		}
		return interpolatedValue;
	}

	private void determineVariogram(VariogramParameters vp, StationsSelection stations) {
		if (!vp.isValid() && inTheoreticalVariogram == null) {
			stations.inNumCloserStations = 0;
			stations.inData = inData;
			VariogramParametersCalculator vpcalCulator = new VariogramParametersCalculator(stations, doDetrended,
					doIncludeZero);
			vpcalCulator.setGlobalVp(vpGlobal);
			vpcalCulator.setGlobalDeTrendedVp(vpGlobalTrend);
			variogramParameters = vpcalCulator.execute();

		} else {
			if (inTheoreticalVariogram != null && inTheoreticalVariogram.get(0)[0] != -9999) {
				variogramParameters = new VariogramParameters.Builder(
						VariogramParameters.getVariogramType(inTheoreticalVariogram.get(5)[0]),
						inTheoreticalVariogram.get(0)[0], inTheoreticalVariogram.get(2)[0],
						inTheoreticalVariogram.get(1)[0]).build();
			} else if (vp.isValid()) {
				variogramParameters = vp;
			} else if (doDetrended && vpGlobalTrend.isValid()) {
				variogramParameters = vpGlobalTrend;
			} else if (vpGlobal.isValid()) {
				variogramParameters = vpGlobal;
			}
		}
	}

	/**
	 * Initializes kriging parameters and verifies input.
	 */
	private VariogramParameters initializeKrigingParameters() {
		VariogramParameters vp = new VariogramParameters.Builder(pSemivariogramType, nugget, range, sill)
				.setLocal(false).setTrend(doDetrended).build();
		if (step == 0) {
			try {
				verifyInput();
				if (!vp.isValid()) {
					this.createDefaulParams();
				}
			} catch (Exception e) {
				pm.errorMessage("Error during input verification: " + e.toString());
			}
		}
		step = step + 1;
		return vp;
	}

	private StationsSelection createAndConfigureStationSelection() {
		StationsSelection stations = new StationsSelection();
		stations.inStations = inStations;
		stations.maxdist = maxdist;
		stations.fStationsid = fStationsid;
		stations.fStationsZ = fStationsZ;
		stations.doLogarithmic = doLogarithmic;
		stations.doIncludezero = doIncludeZero;
		return stations;

	}

	/**
	 * Validates the essential inputs for the kriging model. In detrended mode, both
	 * station and interpolation point elevation fields must be provided.
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
		} else {

		}
	}

	private void createDefaulParams() throws Exception {
		vpGlobal = new VariogramParameters.Builder(globalVariogramType, pNugGlobal, pAGlobal, pSGlobal).build();

		vpGlobalTrend = new VariogramParameters.Builder(globalDetrendedVariogramType, pNugGlobalDeTrended,
				pAGlobalDeTrended, pSGlobalDeTrended).build();
		vpGlobal.setIsLocal(false);
		vpGlobal.setIsTrend(false);
		vpGlobalTrend.setIsLocal(false);
		vpGlobalTrend.setIsTrend(true);
		boolean noEvaluationGlobalParams = (inHValuesPath == null || inHValuesPath.isEmpty());

		if (noEvaluationGlobalParams && !vpGlobal.isValid()) {
			throw new NullPointerException(
					"You provide incomplete fixed parameter for global: nugGlobal, sGlobal, rGlobal or globalVariogramType");
		}

		if (noEvaluationGlobalParams && doDetrended && !vpGlobalTrend.isValid()) {
			throw new NullPointerException(
					"You provide incomplete fixed parameter for global de trended: nugGlobal, sGlobal, rGlobal or globalVariogramType");
		}

		if (!doDetrended) {
			noEvaluationGlobalParams = noEvaluationGlobalParams || vpGlobal.isValid();
		} else {
			noEvaluationGlobalParams = noEvaluationGlobalParams || (vpGlobal.isValid() && vpGlobalTrend.isValid());
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
			vpGlobal = gParam.getGlobalVariogramParameters();
			vpGlobalTrend = gParam.getGlobalVariogramParametersDeTrended();
		}
		if ((nugget == 0 && range == 0 && sill == 0)
				&& (!vpGlobal.isValid() || !(doDetrended && vpGlobalTrend.isValid()))) {
			pm.message("Warning: nugget,range and sill will be evaluate only at each timestep");
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
	 * Calculates the covariance matrix using the theoretical variogram. Each
	 * element [i][j] is computed via TheoreticalVariogram.calculateVGMxyz, which
	 * measures the spatial dissimilarity between points.
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
		outTheoreticalVariogram = variogramParameters.toHashMap();
	}

}
