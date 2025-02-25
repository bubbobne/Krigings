package org.geoframe.blogpost.kriging.rastercase;

import java.awt.image.WritableRaster;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import javax.media.jai.iterator.WritableRandomIter;

import org.geoframe.blogpost.kriging.Kriging;
import org.geoframe.blogpost.kriging.interpolationdata.InterpolationDataProvider;
import org.geoframe.blogpost.kriging.interpolationdata.RasterInterpolationProvider;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.geometry.DirectPosition2D;
import org.hortonmachine.gears.utils.RegionMap;
import org.hortonmachine.gears.utils.coverage.CoverageUtilities;
import org.locationtech.jts.geom.Coordinate;
import org.opengis.geometry.DirectPosition;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import it.geosolutions.jaiext.iterators.RandomIterFactory;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Out;
import oms3.annotations.Status;
/**
 * KrigingRasterCase extends the abstract Kriging class to implement ordinary kriging for raster data.
 * It uses a RasterInterpolationProvider to extract data from a GridCoverage2D input and stores the 
 * interpolated results as a new GridCoverage2D.
 *
 * <p><strong>Key Points and Suggestions:</strong></p>
 * <ul>
 *   <li>
 *     In {@code storeResult()}, the transformation from world coordinates to grid coordinates is performed 
 *     for each interpolated point. If the transformation fails, currently the exception is caught and the stack 
 *     trace is printed. Consider throwing a custom unchecked exception to better manage errors.
 *   </li>
 *   <li>
 *     It is assumed that the ordering of the {@code result} array corresponds exactly to the ordering of 
 *     the keys in {@code interpolatedCoordinatesMap}. Make sure this invariant is maintained.
 *   </li>
 *   <li>
 *     The use of {@code WritableRandomIter} from JAI iterators is appropriate for setting samples in the raster.
 *   </li>
 *   <li>
 *     After updating the raster, a new {@code GridCoverage2D} is built using {@code CoverageUtilities.buildCoverage()}.
 *   </li>
 * </ul>
 */
@Description("Ordinary kriging algorithm.")
@Documentation("Kriging.html")
@Author(name = "Giuseppe Formetta, Daniele Andreis, Silvia Franceschi, Andrea Antonello, Marialaura Bancheri & Francesco Serafin")
@Keywords("Kriging, Hydrology")
@Label("")
@Name("kriging")
@Status()
@License("General Public License Version 3 (GPLv3)")
@SuppressWarnings("nls")
public class KrigingRasterCase extends Kriging {

	@Description("The collection of the points in which the data needs to be interpolated.")
	@In
	public GridCoverage2D inGridCoverage2D = null;

	@Description("The interpolated gridded data ")
	@Out
	public GridCoverage2D outGrid = null;

	/**
    /**
     * Stores the kriging interpolation results in a new GridCoverage2D.
     * The method iterates over the interpolated coordinates, transforms each from world
     * to grid coordinates, and writes the corresponding interpolated value into a writable raster.
     *
     * <p>Potential improvements and error handling:
     * <ul>
     *   <li>
     *     Ensure that the number of interpolated values (in {@code result}) matches the number of points
     *     in {@code interpolatedCoordinatesMap}.
     *   </li>
     *   <li>
     *     Instead of printing the stack trace on transformation errors, consider rethrowing a custom
     *     runtime exception to allow proper error propagation.
     *   </li>
     * </ul>
     * </p>
     *
     * @param result The array of interpolated values.
     * @param interpolatedCoordinatesMap A HashMap mapping unique point IDs to their corresponding Coordinates.
     */
	protected void storeResult(double[] result, HashMap<Integer, Coordinate> interpolatedCoordinatesMap) {
		GridGeometry2D grid = inGridCoverage2D.getGridGeometry();

		final DirectPosition gridPoint = new DirectPosition2D();
		RegionMap regionMap = CoverageUtilities.gridGeometry2RegionParamsMap(grid);
		int cols = regionMap.getCols();
		int rows = regionMap.getRows();

		WritableRaster outWR = CoverageUtilities.createWritableRaster(cols, rows, null, null, null);
		WritableRandomIter outIter = RandomIterFactory.createWritable(outWR, null);

		Set<Integer> pointsToInterpolateIdSett = interpolatedCoordinatesMap.keySet();
		Iterator<Integer> idIterator = pointsToInterpolateIdSett.iterator();
		int c = 0;
		MathTransform transf = grid.getCRSToGrid2D();

		while (idIterator.hasNext()) {
			int id = idIterator.next();
			Coordinate coordinate = (Coordinate) interpolatedCoordinatesMap.get(id);

			DirectPosition point = new DirectPosition2D(grid.getCoordinateReferenceSystem(), coordinate.x,
					coordinate.y);
			try {
				transf.transform(point, gridPoint);
			} catch (MismatchedDimensionException | TransformException e) {
                // SUGGESTION: Instead of only printing the stack trace, consider throwing a custom unchecked exception
                // to allow higher-level error management.
				e.printStackTrace();
			}

			double[] gridCoord = gridPoint.getCoordinate();
			int x = (int) gridCoord[0];
			int y = (int) gridCoord[1];

			outIter.setSample(x, y, 0, result[c]);
			c++;

		}
        // SUGGESTION: Validate that the count 'c' matches the expected number of points to avoid index issues.
		outGrid = CoverageUtilities.buildCoverage("gridded", outWR, regionMap, grid.getCoordinateReferenceSystem());

	}

	@Override
	protected InterpolationDataProvider initializeInterpolatorData() {
		return new RasterInterpolationProvider(inGridCoverage2D);
	}

	

}