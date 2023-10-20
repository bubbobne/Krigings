package krigingsPointCase;

import org.geotools.data.simple.SimpleFeatureCollection;

import experimentalVariogram.ExperimentalVariogram;

public class Utility {

//	public static String[] availableTheorethicalVariogra = new String[] { "exponential", "gaussian", "linear",
//			 "power", "spherical" };
	
	public static String[] availableTheorethicalVariogra = new String[] { "exponential", "linear",
	 "power", "spherical" };
	

	public final static int getVariogramCode(String name) {
		for (int i = 0; i < availableTheorethicalVariogra.length; i++) {
			if (availableTheorethicalVariogra[i] == name) {
				return i;
			}
		}
		return -9999;
	}

	public final static double getLog(double h) {
		return Math.log(h + 1.0);
	}

	public final static double getInverseLog(double h) {
		return Math.exp(h) - 1.0;
	}

	public final static ExperimentalVariogram getExperimentalVariogram(String fStationsid,
			SimpleFeatureCollection inStations, boolean doIncludeZero, int cutoffDivide, double cutoffInput,
			int numCloseStation) {
		ExperimentalVariogram expVariogram = new ExperimentalVariogram();
		expVariogram.fStationsid = fStationsid;
		expVariogram.inStations = inStations;
		expVariogram.doIncludezero = doIncludeZero;
		expVariogram.inNumCloserStations = numCloseStation;
		if (cutoffDivide > 0) {
			expVariogram.Cutoff_divide = cutoffDivide;
		}
		if (cutoffInput > 0) {
			expVariogram.Cutoffinput = cutoffInput;
		}
//		HashMap<Integer, double[]> tmpInData = new HashMap<Integer, double[]>();
//		for (int i = 0; i < idArray.length; i++) {
//			tmpInData.put(idArray[i], new double[] { hresiduals[i] });
//		}
//		expVariogram.inData = tmpInData;
		return expVariogram;
	}
}
