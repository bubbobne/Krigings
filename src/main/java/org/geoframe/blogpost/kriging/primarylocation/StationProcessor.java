package org.geoframe.blogpost.kriging.primarylocation;

import java.util.Arrays;
import java.util.HashMap;

import org.geoframe.blogpost.kriging.pointcase.Kriging;
import org.locationtech.jts.geom.Coordinate;


public class StationProcessor {
	private StationsSelection stations;
	private double[] xStations;
	private double[] yStations;
	private double[] zStations;
	private double[] hStations;
	private double[] hResiduals;
	private int count;
	private double trendCoeff;
	private double trendIntercept;
	private boolean areAllEquals;

	public StationProcessor(StationsSelection stations) {
		this.stations = stations;
	}
	
	/**
	 * Execute station selection for the given coordinate.
	 */
	public void updateForCoordinate(Coordinate coordinate, HashMap<Integer, double[]> inData, int inNumCloserStations,
			double maxdist, boolean doDetrended) throws Exception {
		if (inNumCloserStations > 0) {
			stations.inNumCloserStations = inNumCloserStations;
		}
		if (maxdist > 0) {
			stations.maxdist = maxdist;
		}
		if (coordinate!=null && (maxdist > 0 || inNumCloserStations > 0)) {
			stations.idx = coordinate.x;
			stations.idy = coordinate.y;
		}
		stations.inData = inData;
		stations.execute();

		this.xStations = stations.xStationInitialSet;
		this.yStations = stations.yStationInitialSet;
		this.zStations = stations.zStationInitialSet;
		this.hStations = stations.hStationInitialSet;
		this.count = xStations.length - 1;
		this.areAllEquals = stations.areAllEquals;

		// Compute residuals and trends if applicable
		if (this.count > 0) {
			ResidualsEvaluator residualsEvaluator = getResidualsEvaluator(
					Arrays.copyOfRange(this.zStations, 0, this.count),
					Arrays.copyOfRange(this.hStations, 0, this.count), doDetrended);
			this.hResiduals = residualsEvaluator.hResiduals;
			this.trendCoeff = residualsEvaluator.trend_coefficient;
			this.trendIntercept = residualsEvaluator.trend_intercept;
		}
	}

	// Getters for the various station data fields
	public double[] getXStations() {
		return xStations;
	}

	public double[] getYStations() {
		return yStations;
	}

	public double[] getZStations() {
		return zStations;
	}

	public double[] getHStations() {
		return hStations;
	}

	public double[] getHResiduals() {
		return hResiduals;
	}

	public int getCount() {
		return count;
	}

	public double getTrendCoeff() {
		return trendCoeff;
	}

	public double getTrendIntercept() {
		return trendIntercept;
	}

	public boolean areAllEquals() {
		return areAllEquals;
	}

	private ResidualsEvaluator getResidualsEvaluator(double[] zStations, double[] hStations, boolean doDetrended) {
		ResidualsEvaluator residualsEvaluator = new ResidualsEvaluator();
		residualsEvaluator.doDetrended = doDetrended;
		residualsEvaluator.hStations = hStations;
		residualsEvaluator.zStations = zStations;
		residualsEvaluator.regressionOrder = Kriging.REGRESSION_ORDER;
		residualsEvaluator.process();
		return residualsEvaluator;
	}

}
