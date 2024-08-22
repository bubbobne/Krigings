package theoreticalVariogram;

import java.util.HashMap;

import org.hortonmachine.gears.libs.modules.HMConstants;

import krigingsPointCase.Utility;

public class VariogramParameters {
	private double nugget;
	private double range;
	private double sill;
	private String modelName =null;
	private boolean isTrend;

	public boolean getIsTrend() {
		return isTrend;
	}

	public void setIsTrend(boolean isTrend) {
		this.isTrend = isTrend;
	}

	private boolean isLocal;

	public boolean getIsLocal() {
		return isLocal;
	}

	public void setIsLocal(boolean isLocal) {
		this.isLocal = isLocal;
	}

	public VariogramParameters(String modelName, double nugget, double range, double sill) {
		this.modelName = modelName;
		this.nugget = nugget;
		this.range = range;
		this.sill = sill;
	}

	public VariogramParameters() {
		// TODO Auto-generated constructor stub
	}

	public double getNugget() {
		return nugget;
	}

	public double getRange() {
		return range;
	}

	public double getSill() {
		return sill;
	}

	public String getModelName() {
		return modelName;
	}

	public HashMap<Integer, double[]> toHashMap() {
		HashMap<Integer, double[]> outVariogramParams = new HashMap<Integer, double[]>();
		if (modelName != null) {
			outVariogramParams.put(0, new double[] { nugget });
			outVariogramParams.put(1, new double[] { sill });
			outVariogramParams.put(2, new double[] { range });
			outVariogramParams.put(3, new double[] { isLocal ? 0.0 :1.0 });
			outVariogramParams.put(4, new double[] { isTrend ? 0.0 :1.0 });
			outVariogramParams.put(5, new double[] { Utility.getVariogramCode(modelName) });
		} else {
			outVariogramParams = new HashMap<Integer, double[]>();
			outVariogramParams.put(0, new double[] { HMConstants.doubleNovalue });
			outVariogramParams.put(1, new double[] { HMConstants.doubleNovalue });
			outVariogramParams.put(2, new double[] { HMConstants.doubleNovalue });
			outVariogramParams.put(3, new double[] { HMConstants.doubleNovalue });
			outVariogramParams.put(4, new double[] { HMConstants.doubleNovalue });
			outVariogramParams.put(4, new double[] { HMConstants.doubleNovalue });
		}
		return outVariogramParams;
	}

}