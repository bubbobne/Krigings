package KrigingsTests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.DataUtilities;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.filter.text.cql2.CQLException;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.CRS;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.hortonmachine.gears.libs.monitor.DummyProgressMonitor;
import org.hortonmachine.gears.utils.coverage.CoverageUtilities;
import org.hortonmachine.gears.utils.features.FeatureUtilities;
//import org.hortonmachine.hmachine.utils.HMTestMaps;
import org.junit.Before;
import org.junit.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.geometry.MismatchedDimensionException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.TransformException;

import krigingsPointCase.Krigings;

/**
 * Test the Krigings model.
 * 
 * @author daniele andreis
 * 
 */
public class TestKriging {

	private File stazioniFile;
	private File puntiFile;
	private File KrigingsRainFile;
	private String interpolatedRainPath;
	private File KrigingsRain2File;
	private File KrigingsRain3File;
	private File krigingRain4File;
	private File stazioniGridFile;

	@Before
	public void setUp() throws Exception {

		URL stazioniUrl = this.getClass().getClassLoader().getResource("rainstations.shp");
		stazioniFile = new File(stazioniUrl.toURI());

		URL puntiUrl = this.getClass().getClassLoader().getResource("basins_passirio_width0.shp");
		puntiFile = new File(puntiUrl.toURI());

		URL KrigingsRainUrl = this.getClass().getClassLoader().getResource("rain_test.csv");
		KrigingsRainFile = new File(KrigingsRainUrl.toURI());

		URL KrigingsRain2Url = this.getClass().getClassLoader().getResource("rain_test2A.csv");
		KrigingsRain2File = new File(KrigingsRain2Url.toURI());

		URL KrigingsRain3Url = this.getClass().getClassLoader().getResource("rain_test3A.csv");
		KrigingsRain3File = new File(KrigingsRain3Url.toURI());

		URL stazioniGridUrl = this.getClass().getClassLoader().getResource("rainstationgrid.shp");
		stazioniGridFile = new File(stazioniGridUrl.toURI());

		URL KrigingsRain4Url = this.getClass().getClassLoader().getResource("rain_test_grid.csv");
		krigingRain4File = new File(KrigingsRain4Url.toURI());

		File interpolatedRainFile = new File(KrigingsRainFile.getParentFile(), "Krigings_interpolated.csv");
		interpolatedRainPath = interpolatedRainFile.getAbsolutePath();
		// interpolatedRainPath = interpolatedRainPath.replaceFirst("target",
		// "src" + File.separator + File.separator + "test");
		interpolatedRainPath = interpolatedRainPath.replaceFirst("target", "src" + File.separator + "test");

		interpolatedRainPath = interpolatedRainPath.replaceFirst("test-classes", "resources");

	}

	@Test
	public void testKriging() throws URISyntaxException, SchemaException, CQLException {
		//
		String stationIdField = "ID_PUNTI_M";
		System.out.println("Working Directory = " + System.getProperty("user.dir"));

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		try {
			stationsReader.readFeatureCollection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		// OmsShapefileFeatureReader interpolatedPointsReader = new
		// OmsShapefileFeatureReader();
		// interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		// interpolatedPointsReader.readFeatureCollection();

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = krigingRain4File.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		reader.initProcess();

		String inPathToDistance = "src/test/resources/distance.csv";
		String inPathToVariance = "src/test/resources/exp_var.csv";

		int timeStepMinutes = 60;
		String fId = "ID";




		Krigings kriging = new Krigings();
		GeometryFactory geomFactory = new GeometryFactory();
		Point point = (Point) geomFactory.createPoint(new Coordinate(0.0, 50.0));
		SimpleFeatureType TYPE =
	                DataUtilities.createType(
	                        "Location",
	                        "the_geom:Point:srid=32632,"
	                                + // <- the geometry attribute: Point type
	                                "ID:int,"

	                        );
		
		
		Filter filter = CQL.toFilter(stationIdField + " = 1331");
		SimpleFeatureCollection subCollection = stationsFC.subCollection(filter);
		assertTrue(subCollection.size() == 1);

		SimpleFeature station = subCollection.features().next();
		Geometry geometry = (Geometry) station.getDefaultGeometry();
		Coordinate stationCoordinate = geometry.getCoordinate();
        SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);
        featureBuilder.add(point);
        featureBuilder.add(1);
        
        SimpleFeature feature = featureBuilder.buildFeature(null);

        List<SimpleFeature> features = new ArrayList<>();
        features.add(feature);
        SimpleFeatureCollection collection = new ListFeatureCollection(TYPE, features);
		kriging.inInterpolate = subCollection;
		kriging.pSemivariogramType = "exponential";

		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;

		// kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = stationIdField;

		// it doesn't execute the model with log value.
		/*
		 * Set up the model in order to use the variogram with an explicit integral
		 * scale and variance.
		 */
		// kriging.pVariance = 3.5;
		// kriging.pIntegralscale = new double[]{10000, 10000, 100};

//        kriging.pA = 123537.0;
//        kriging.pNug = 0.0;
//        kriging.pS = 1.678383;
		/*
		 * Set up the model in order to run with a FeatureCollection as point to
		 * interpolated. In this case only 2D.
		 */

		// OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		// writer.file = interpolatedRainPath;
		//
		// writer.tStart = reader.tStart;
		// writer.tTimestep = reader.tTimestep;
		while (reader.doProcess) {
			try {
				reader.nextRecord();

				HashMap<Integer, double[]> id2ValueMap = reader.outData;
				kriging.inData = id2ValueMap;


				kriging.executeKriging();
				/*
				 * Extract the result.
				 */

				double[] values = id2ValueMap.get(1331);
			

				HashMap<Integer, double[]> result = kriging.outData;
				new Point2D.Double(stationCoordinate.x, stationCoordinate.y);
				Set<Integer> pointsToInterpolateResult = result.keySet();
				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
				int id = iteratorTest.next();
				double[] actual = result.get(id);
				assertEquals(actual[0], values[0], 0.01);

				// HashMap<Integer, double[]> result = kriging.outFolder;
				// Set<Integer> pointsToInterpolateResult = result.keySet();
				// Iterator<Integer> iteratorTest = pointsToInterpolateResult
				// .iterator();

				int iii = 0;
				// while (iteratorTest.hasNext() && iii<12) {
				// double expected;
				// if (j == 0) {
				// expected = 0.3390869;
				// } else if (j == 1) {
				// expected = 0.2556174;
				// } else if (j == 2) {
				// expected = 0.2428944;
				// } else if (j == 3) {
				// expected = 0.2613782;
				// } else if (j == 4) {
				// expected = 0.3112850;
				// } else if (j == 5) {
				// expected = 0.2983679;
				// } else if (j == 6) {
				// expected = 0.3470377;
				// } else if (j == 7) {
				// expected = 0.3874065;
				// } else if (j == 8) {
				// expected = 0.2820323;
				// } else if (j == 9) {
				// expected = 0.1945515;
				// } else if (j == 10) {
				// expected = 0.1698022;
				// } else if (j == 11) {
				// expected = 0.2405134;
				// } else if (j == 12) {
				// expected = 0.2829313;
				// } else {
				// expected = 1.0;
				// }				// expected = 0.1698022;

				//

				//
				// int id = iteratorTest.next();
				// double[] actual = result.get(id);
				// iii+=1;
				//
				// //assertEquals(expected, actual[0], 0.001);
				// j=j+1;
				// }
				// iii=0;
				// j=0;
				// writer.inData = result;
				// writer.writeNextLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (CQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	// writer.close();
	}

	//
	//

	// ///////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////////TEST 1
	// PASSA////////////////////////////////////////////////////
	// /////////////////////////////////////////////////////////////////////////////////////////
	@Test
	public void testKrigings1() throws Exception {

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = KrigingsRainFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";

		reader.initProcess();

		Krigings kriging = new Krigings();
		kriging.pm = new DummyProgressMonitor();

		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";

		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";
		kriging.cutoffDivide = 10;
		kriging.inHValuesPath = KrigingsRainFile.getAbsolutePath();
		// it d		kriging.pSemivariogramType = "exponential";
//oesn't execute the model with log value.
		/*
		 * Set up the model in order to use the variogram with an explicit integral
		 * scale and variance.
		 */
		// Krigings.pVariance = 3.5;
		// Krigings.pIntegralscale = new double[]{10000, 10000, 100};
//		Krigings.range = 123537.0;
//		Krigings.nugget = 0.0;
//		Krigings.sill = 1.678383;
		/*
		 * Set up the model in order to run with a FeatureCollection as point to
		 * interpolated. In this case only 2D.
		 */
		kriging.pSemivariogramType = "exponential";

		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = interpolatedRainPath;

		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		int j = 0;
		while (reader.doProcess) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = kriging.outData;
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();

			int iii = 0;

			while (iteratorTest.hasNext() && iii < 12) {
				double expected;
				if (j == 0) {
					expected = 0.3390869;
				} else if (j == 1) {
					expected = 0.2556174;
				} else if (j == 2) {
					expected = 0.2428944;
				} else if (j == 3) {
					expected = 0.2613782;
				} else if (j == 4) {
					expected = 0.3112850;
				} else if (j == 5) {
					expected = 0.2983679;
				} else if (j == 6) {
					expected = 0.3470377;
				} else if (j == 7) {
					expected = 0.3874065;
				} else if (j == 8) {
					expected = 0.2820323;
				} else if (j == 9) {
					expected = 0.1945515;
				} else if (j == 10) {
					expected = 0.1698022;
				} else if (j == 11) {
					expected = 0.2405134;
				} else if (j == 12) {
					expected = 0.2829313;
				} else {
					expected = 1.0;
				}

				int id = iteratorTest.next();
				double[] actual = result.get(id);
				iii += 1;

				assertEquals(expected, actual[0], 0.001);
				System.out.println(" " + expected + "  " + actual[0]);
				j = j + 1;
			}
			iii = 0;
			j = 0;
			writer.inData = result;
			writer.writeNextLine();

		}

		reader.close();
		writer.close();
	}

	// ///////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////////FINE TEST
	// 1PASSA////////////////////////////////////////////////////
	// /////////////////////////////////////////////////////////////////////////////////////////
	//
	//
	//
	// ///////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////////TEST 2
	// PASSA////////////////////////////////////////////////////
	// /////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * Run the Krigings models.
	 *
	 * <p>
	 * This is the case which all the station have the same value.
	 * </p>
	 * 
	 * @throws Exception
	 * @throws Exception
	 */
	@Test
	public void testKrigings2() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = KrigingsRain2File.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";

		reader.initProcess();

		Krigings Krigings = new Krigings();
		Krigings.pm = new DummyProgressMonitor();

		Krigings.inStations = stationsFC;
		Krigings.fStationsid = "ID_PUNTI_M";

		Krigings.inInterpolate = interpolatedPointsFC;
		Krigings.fInterpolateid = "netnum";

		// it doesn't execute the model with log value.
		/*
		 * Set up the model in order to use the variogram with an explicit integral
		 * scale and variance.
		 */

		/*
		 * Set up the model in order to run with a FeatureCollection as point to
		 * interpolated. In this case only 2D.
		 */

		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = interpolatedRainPath;

		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;

		while (reader.doProcess) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			Krigings.inData = id2ValueMap;
			Krigings.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = Krigings.outData;
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iterator = pointsToInterpolateResult.iterator();
			while (iterator.hasNext()) {
				int id = iterator.next();
				double[] actual = result.get(id);
				assertEquals(1.0, actual[0], 0);
			}
			writer.inData = result;
			writer.writeNextLine();
		}

		reader.close();
		writer.close();
	}
	// /////////////////////////////////////////////////////////////////////////////////////////
	// ///////////////////////////////FINE TEST 2
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////

	// /////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////// TEST 3
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	// /**
	// * Run the Krigings models.
	// *
	// * <p>
	// * This is the case that defaultMode=0.
	// * </p>
	// * @throws Exception
	// * @throws Exception
	// */
	@Test
	public void testKrigings4() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = KrigingsRainFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";

		reader.initProcess();

		Krigings kriging = new Krigings();
		kriging.pm = new DummyProgressMonitor();

		kriging.inStations = stationsFC;
		kriging.fStationsid = "ID_PUNTI_M";

		kriging.inInterpolate = interpolatedPointsFC;
		kriging.fInterpolateid = "netnum";

		// it doesn't execute the model with log value.
		/*
		 * Set up the model in order to use the variogram with an explicit integral
		 * scale and variance.
		 */
		/*
		 * Set up the model in order to run with a FeatureCollection as point to
		 * interpolated. In this case only 2D.
		 */
		kriging.pSemivariogramType = "exponential";
		kriging.cutoffDivide = 3;
		kriging.doIncludeZero = false;
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = interpolatedRainPath;

		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;

		while (reader.doProcess) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = kriging.outData;
//			double[][] test = HMTestMaps.outKriging4;
//			for (int i = 0; i < test.length; i++) {
//				double actual = result.get((int) test[i][0])[0];
//				double expected = test[i][1];
//				assertEquals(expected, actual, 0.01);
//			}

			writer.inData = result;
			writer.writeNextLine();
		}

		reader.close();
		writer.close();
	}
	// /////////////////////////////////////////////////////////////////////////////////////////
	// ///////////////////////////////FINE TEST 3
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////

	// /////////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////// TEST 4
	// PASSA////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Run the Krigings models.
	 *
	 * <p>
	 * This is the case which there is only one station.
	 * </p>
	 * 
	 * @throws Exception
	 * @throws Exception
	 */
	@Test
	public void testKrigings5() throws Exception {
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		interpolatedPointsReader.readFeatureCollection();
		SimpleFeatureCollection interpolatedPointsFC = interpolatedPointsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = KrigingsRain3File.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";

		reader.initProcess();

		Krigings Krigings = new Krigings();
		Krigings.pm = new DummyProgressMonitor();

		Krigings.inStations = stationsFC;
		Krigings.fStationsid = "ID_PUNTI_M";

		Krigings.inInterpolate = interpolatedPointsFC;
		Krigings.fInterpolateid = "netnum";

		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = interpolatedRainPath;

		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		int j = 0;
		while (reader.doProcess) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			Krigings.inData = id2ValueMap;
			Krigings.executeKriging();
			/*
			 * Extract the result.
			 */
			HashMap<Integer, double[]> result = Krigings.outData;
			Set<Integer> pointsToInterpolateResult = result.keySet();
			Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
			double expected;
			if (j == 0) {
				expected = 10.0;
			} else if (j == 1) {
				expected = 15;
			} else if (j == 2) {
				expected = 1;
			} else if (j == 3) {
				expected = 2;
			} else if (j == 4) {
				expected = 2;
			} else if (j == 5) {
				expected = 0;
			} else if (j == 6) {
				expected = 0;
			} else if (j == 7) {
				expected = 23;
			} else if (j == 8) {
				expected = 50;
			} else if (j == 9) {
				expected = 70;
			} else if (j == 10) {
				expected = 30;
			} else if (j == 11) {
				expected = 10;
			} else if (j == 12) {
				expected = 2;
			} else {
				expected = 1.0;
			}

			while (iteratorTest.hasNext()) {
				int id = iteratorTest.next();
				double[] actual = result.get(id);

				assertEquals(expected, actual[0], 0);
			}
			writer.inData = result;
			writer.writeNextLine();
			j++;
		}

		reader.close();
		writer.close();
	}

	// ///////////////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////FINE TEST 4
	// PASSA////////////////////////////////////////////////////
	// /////////////////////////////////////////////////////////////////////////////////////

	private OmsTimeSeriesIteratorReader getTimeseriesReader(String inPath, String id, String startDate,
			int timeStepMinutes) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
