package theoreticalVariogram;

import java.util.Arrays;
import java.util.HashMap;

import experimentalVariogram.ExperimentalVariogram;
import krigingsPointCase.Krigings;
import krigingsPointCase.ResidualsEvaluator;
import krigingsPointCase.StationsSelection;

public class VariogramParametersCalculator {
	private VariogramParameters globalVP = new VariogramParameters();
	private VariogramParameters globalDeTrendedVP = new VariogramParameters();
	private StationsSelection stations = null;
	private boolean doDetrend = false;
	private String semivariogramType = null;
	private int cutoffDivide;
	private double cutoffInput;

	public VariogramParametersCalculator(StationsSelection stations, boolean doDetrend, Boolean doIncludeZero) {
		this.stations = stations;
		this.doDetrend = doDetrend;
		try {
			int tempN = 0;
			if (stations.inNumCloserStations != 0) {
				tempN = stations.inNumCloserStations;
				stations.inNumCloserStations = 0;
			}
			this.stations.execute();
			stations.inNumCloserStations = tempN;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public VariogramParameters execute() {
		if (!stations.areAllEquals && stations.n1 > 2) {
			int n1 = stations.hStationInitialSet.length - 1;
			n1 = stations.hStationInitialSet.length - 1;
			ResidualsEvaluator rEvaluator = new ResidualsEvaluator();
			rEvaluator.doDetrended = doDetrend;
			rEvaluator.hStations = Arrays.copyOfRange(stations.hStationInitialSet, 0, n1);
			rEvaluator.zStations = Arrays.copyOfRange(stations.zStationInitialSet, 0, n1);
			rEvaluator.regressionOrder = Krigings.REGRESSION_ORDER;
			rEvaluator.process();
			int[] idStations = Arrays.copyOfRange(stations.idStationInitialSet, 0, n1);
			double[] hResiduals = rEvaluator.hResiduals;
			hResiduals = rEvaluator.hResiduals;
			VariogramParamsEvaluator variogramParamsEvaluator = new VariogramParamsEvaluator();
			variogramParamsEvaluator.expVar = getExperimentalVariogram(hResiduals, idStations);
			variogramParamsEvaluator.pSemivariogramType = this.semivariogramType;
			variogramParamsEvaluator.proces();
			boolean variogramOk = variogramParamsEvaluator.nugget >= 0 && variogramParamsEvaluator.sill > 0
					&& variogramParamsEvaluator.range > 0 && variogramParamsEvaluator.isFitGood;
			if (variogramOk || (!globalVP.isValid()) || (doDetrend && globalDeTrendedVP.isValid())) {
				VariogramParameters myVariogramParam = new VariogramParameters(
						variogramParamsEvaluator.outSemivariogramType, variogramParamsEvaluator.nugget,
						variogramParamsEvaluator.range, variogramParamsEvaluator.sill);
				myVariogramParam.setIsLocal(true);
				myVariogramParam.setIsTrend(rEvaluator.isPValueOk);
				return myVariogramParam;
			}
		}
		if (doDetrend && globalDeTrendedVP.isValid()) {
			return globalDeTrendedVP;
		} else if (globalVP.isValid()) {
			return globalVP;
		}
		return null;
	}

	private ExperimentalVariogram getExperimentalVariogram(double[] hresiduals, int[] idArray) {
		ExperimentalVariogram expVariogram = ExperimentalVariogram.create(stations.fStationsid, stations.inStations,
				stations.doIncludezero, cutoffDivide, cutoffInput, 0);
		HashMap<Integer, double[]> tmpInData = new HashMap<Integer, double[]>();
		for (int i = 0; i < idArray.length; i++) {
			tmpInData.put(idArray[i], new double[] { hresiduals[i] });
		}
		expVariogram.inData = tmpInData;
		return expVariogram;
	}

	public void setGlobalDeTrendedVp(VariogramParameters vpGlobalDetrended) {
		// TODO Auto-generated method stub
		if (vpGlobalDetrended != null)
			this.globalDeTrendedVP = vpGlobalDetrended;
	}

	public void setGlobalVp(VariogramParameters vpGlobal) {
		// TODO Auto-generated method stub
		if (vpGlobal != null)
			this.globalVP = vpGlobal;
	}
}
