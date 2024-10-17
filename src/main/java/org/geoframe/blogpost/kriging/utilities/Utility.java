package org.geoframe.blogpost.kriging.utilities;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.hortonmachine.gears.libs.monitor.IHMProgressMonitor;
import org.hortonmachine.hmachine.i18n.HortonMessageHandler;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;

public class Utility {

//	public static String[] availableTheorethicalVariogra = new String[] { "exponential", "gaussian", "linear",
//			 "power", "spherical" };



	public final static double getLog(double h) {
		return Math.log(h + 1.0);
	}

	public final static HashMap<Integer, double[]> getLog(HashMap<Integer, double[]> h) {
		HashMap<Integer, double[]> logMap = new HashMap<>();
		for (Map.Entry<Integer, double[]> entry : h.entrySet()) {
			double[] originalArray = entry.getValue();
			double[] logArray = new double[originalArray.length];
			for (int i = 0; i < originalArray.length; i++) {
				if (originalArray[i] >= 0) {
					logArray[i] = getLog(originalArray[i]);
				} else {
					// Handle non-positive values, e.g., by setting to NaN or skipping
					logArray[i] = Double.NaN; // Or you can use Double.NEGATIVE_INFINITY
				}
			}
			logMap.put(entry.getKey(), logArray);
		}

		// Optionally, you can now replace the old map with the new map if modification
		// in place isn't desired
		return logMap;
	}

	public final static double getInverseLog(double h) {
		return Math.exp(h) - 1.0;
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
	public final static LinkedHashMap<Integer, Coordinate> getCoordinate(SimpleFeatureCollection collection,
			String idField, String fPointZ, IHMProgressMonitor pm, HortonMessageHandler msg) throws Exception {
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

}
