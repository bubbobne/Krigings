package KrigingsTests;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geoframe.blogpost.kriging.pointcase.KrigingPointCase;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.junit.Test;
import org.opengis.filter.Filter;

public class TestValDiNon {
	@Test
	public void testKrigingSic97WithTrend() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "ID";
		// 100 station to training model
	    URL stationShpUrl = this.getClass().getClassLoader().getResource("MeteoStations.shp");
		File stazioniGridFile = new File(stationShpUrl.getFile());
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		URL observedRain4Url = this.getClass().getClassLoader().getResource("precipitation_cleaned.csv");
		File observedFile = new File(observedRain4Url.toURI());
		reader.file = observedFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2017-04-12 10:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2017-12-12 05:00";
		reader.fileNovalue = "-9999.0";
		reader.initProcess();

		String fId = "ID";
		KrigingPointCase kriging = new KrigingPointCase();

		int idToCheck = 25;
		Filter filter = CQL.toFilter("ID = " + idToCheck);

		kriging.inInterpolate = stationsFC.subCollection(filter);
		//kriging.pSemivariogramType = "exponential";
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
		kriging.fInterpolateid = stationIdField;
//        kriging.inHValuesPath= observedFile.getAbsolutePath();
//		kriging.cutoffDivide = 10;
		kriging.doDetrended = true;
		kriging.fPointZ = "z_dem";
		kriging.fStationsZ = "z_dem";
//		kriging.inNumCloserStations = 5;
		kriging.doLogarithmic = false;
		kriging.boundedToZero = true;

//		kriging.globalVariogramType  ="exponential";
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/krigings/PointCase/stations_search_" + idToCheck + ".csv";
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
		OmsTimeSeriesIteratorWriter writer2 = new OmsTimeSeriesIteratorWriter();
		writer2.file = "resources/Output/krigings/PointCase/params_" + idToCheck + ".csv";
		writer2.tStart = reader.tStart;
		writer2.tTimestep = reader.tTimestep;
		
		
		while (reader.doProcess) {
			try {
				reader.nextRecord();
				System.out.println("new T:" + reader.tCurrent);
				HashMap<Integer, double[]> id2ValueMap = reader.outData;

				HashMap<Integer, double[]> predictedGstatR = new HashMap<Integer, double[]>();
				predictedGstatR.put(idToCheck, id2ValueMap.get(idToCheck));
				id2ValueMap.remove(idToCheck);
				kriging.inData = id2ValueMap;
				kriging.executeKriging();
				HashMap<Integer, double[]> result = kriging.outData;
				Set<Integer> pointsToInterpolateResult = result.keySet();
				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
				while (iteratorTest.hasNext()) {
					int id = iteratorTest.next();
					double[] values = result.get(id);
					double[] actual = predictedGstatR.get(id);
					assertEquals(actual[0], values[0], 2);
					System.out.println("actual is:" + actual[0] + " evaluate " + values[0]);
				}

				writer.inData = kriging.outData;
//				writer2.inData = kriging.outVariogramParams;
				writer2.writeNextLine();

				writer.writeNextLine();
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
		writer.close();
		writer2.close();

		reader.close();

		// writer.close();
	}

	@Test
	public void testKrigingBasin() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "ID";
		// 100 station to training model
		URL stationShpUrl = this.getClass().getClassLoader().getResource("MeteoStations.shp");
		File stazioniGridFile = new File(stationShpUrl.toURI());
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		URL observedRain4Url = this.getClass().getClassLoader().getResource("precipitation_cleaned.csv");
		File observedFile = new File(observedRain4Url.toURI());
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = observedFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2022-09-14 19:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2017-12-12 05:00";
		reader.fileNovalue = "-9999.0";
		reader.initProcess();

		String fId = "ID";
		KrigingPointCase kriging = new KrigingPointCase();

		URL predictUrl = this.getClass().getClassLoader().getResource("centroid_ID_1881.shp");
		File predictFile = new File(predictUrl.toURI());
		OmsShapefileFeatureReader predictReader = new OmsShapefileFeatureReader();
		predictReader.file = predictFile.getAbsolutePath();
		predictReader.readFeatureCollection();
		SimpleFeatureCollection predictFC = predictReader.geodata;
		kriging.inInterpolate = predictFC;
		kriging.pSemivariogramType = "exponential";
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
		kriging.fInterpolateid = "basinid";
//		kriging.inHValuesPath = observedFile.getAbsolutePath();
//		kriging.cutoffDivide = 15;
		kriging.doDetrended = true;
		kriging.fPointZ = "elev_m";
		kriging.fStationsZ = "z_dem";
//		kriging.inNumCloserStations = 5;
		kriging.doLogarithmic = false;
		kriging.boundedToZero = true;

//		kriging.inHValuesPath = observedFile.getAbsolutePath();
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.tStart = reader.tStart;
		writer.file = "resources/Output/krigings/PointCase/basin1881_2.csv";

		writer.tTimestep = reader.tTimestep;
		
		

		while (reader.doProcess) {
			try {
				reader.nextRecord();
				System.out.println("new T:" + reader.tCurrent);
				HashMap<Integer, double[]> id2ValueMap = reader.outData;
				;
				kriging.inData = id2ValueMap;
				kriging.executeKriging();
				HashMap<Integer, double[]> result = kriging.outData;
				Set<Integer> pointsToInterpolateResult = result.keySet();
				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
				while (iteratorTest.hasNext()) {
					int id = iteratorTest.next();
					double[] values = result.get(id);
					System.out.println(values[0]);
					//assertEquals(actual[0], values[0], 2);
				}

				writer.inData = kriging.outData;
				writer.writeNextLine();
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
		writer.close();

		reader.close();

		// writer.close();
	}

//	@Test
//	public void testKrigingOmsConsole() throws URISyntaxException, SchemaException, CQLException, IOException {
//		//
//		String stationIdField = "ID";
//		// 100 station to training model
//		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
//		String home = "/home/andreisd/Documents/project/GWS2022/OMS_project_val_di_non/OMS_project";
//		stationsReader.file = home + "/data/Trentino/Noce/MeteoStations.shp";
//		stationsReader.readFeatureCollection();
//		SimpleFeatureCollection stationsFC = stationsReader.geodata;
//
//		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
//		reader.file = home + "/data/Trentino/Noce/precipitation_normal_score.csv";
//		reader.idfield = "ID";
//		reader.tStart = "2017-11-30 01:00";
//		reader.tTimestep = 60;
//		reader.tEnd = "2017-12-01 05:00";
//		reader.fileNovalue = "-9999.0";
//		reader.initProcess();
//
//		String fId = "ID";
//		Krigings kriging = new Krigings();
//
//		int idToCheck = 25;
//		OmsShapefileFeatureReader predictReader = new OmsShapefileFeatureReader();
//		predictReader.file = home
//				+ "/output/geomorphology/Noce_network_50/km10/geoframe/small_basin/1491006/centroid_ID_1491006.shp";
//		predictReader.readFeatureCollection();
//		SimpleFeatureCollection predictFC = predictReader.geodata;
//
//		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
//		writer.file = "resources/Output/krigings/PointCase/stations_" + idToCheck + ".csv";
//		writer.tStart = reader.tStart;
//		writer.tTimestep = reader.tTimestep;
//		while (reader.doProcess) {
//			try {
//
//				kriging.inInterpolate = predictFC;
//				kriging.pSemivariogramType = "exponential";
//				kriging.inStations = stationsFC;
//				kriging.fStationsid = stationIdField;
//				kriging.fInterpolateid = "basinid";
//				kriging.inHValuesPath = home + "/data/Trentino/Noce/precipitation_normal_score.csv";
//				kriging.cutoffDivide = 10;
//				kriging.doDetrended = false;
//				kriging.fPointZ = "elev_m";
//				kriging.fStationsZ = "z_dem";
////				kriging.inNumCloserStations = 5;
//				kriging.doLogarithmic = false;
//				kriging.boundedToZero = false;
//				reader.nextRecord();
//				System.out.println("new T:" + reader.tCurrent);
//				HashMap<Integer, double[]> id2ValueMap = reader.outData;
//
//				HashMap<Integer, double[]> predictedGstatR = new HashMap<Integer, double[]>();
//				predictedGstatR.put(idToCheck, id2ValueMap.get(idToCheck));
//				id2ValueMap.remove(idToCheck);
//				kriging.inData = id2ValueMap;
//				kriging.executeKriging();
//				HashMap<Integer, double[]> result = kriging.outData;
//				Set<Integer> pointsToInterpolateResult = result.keySet();
//				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
//				while (iteratorTest.hasNext()) {
//					int id = iteratorTest.next();
//					double[] values = result.get(id);
//					// assertEquals(actual[0], values[0], 2);
//				}
//
//				writer.inData = kriging.outData;
//				writer.writeNextLine();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			} catch (CQLException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			} catch (Exception e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		writer.close();
//
//		reader.close();
//
//		// writer.close();
//	}

}
