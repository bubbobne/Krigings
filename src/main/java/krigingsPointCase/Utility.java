package krigingsPointCase;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.geotools.data.simple.SimpleFeatureCollection;

import experimentalVariogram.ExperimentalVariogram;

public class Utility {

//	public static String[] availableTheorethicalVariogra = new String[] { "exponential", "gaussian", "linear",
//			 "power", "spherical" };
	
	public static String[] availableTheorethicalVariogra = new String[] { "exponential", "linear", "power", "spherical"};
	

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

	public final static HashMap<Integer, double[]> getLog(HashMap<Integer, double[]> h ) { 
		HashMap<Integer, double[]> logMap = new HashMap<>();
        for (Map.Entry<Integer, double[]> entry : h.entrySet()) {
            double[] originalArray = entry.getValue();
            double[] logArray = new double[originalArray.length];
            for (int i = 0; i < originalArray.length; i++) {
                if (originalArray[i] >= 0) {
                    logArray[i] = getLog(originalArray[i]);
                } else {
                    // Handle non-positive values, e.g., by setting to NaN or skipping
                    logArray[i] = Double.NaN;  // Or you can use Double.NEGATIVE_INFINITY
                }
            }
            logMap.put(entry.getKey(), logArray);
        }

        // Optionally, you can now replace the old map with the new map if modification in place isn't desired
        return logMap;
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
