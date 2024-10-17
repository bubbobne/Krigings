package KrigingsTests;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geoframe.blogpost.kriging.pointcase.Kriging;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.junit.Test;

public class TestKrigingStandardFromR {

	@Test
	public void testKrigingSic97() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "id";
		// 100 station to training model
		URL stazioniGridUrl = this.getClass().getClassLoader().getResource("observed.shp");
		File stazioniGridFile = new File(stazioniGridUrl.toURI());
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		URL observedRain4Url = this.getClass().getClassLoader().getResource("observed_H.csv");
		File observedFile = new File(observedRain4Url.toURI());
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = observedFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2022-12-06 17:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		reader.initProcess();
		OmsTimeSeriesIteratorReader testReaderValue = new OmsTimeSeriesIteratorReader();
		URL testRain4Url = this.getClass().getClassLoader().getResource("test_H.csv");
		File testFile = new File(testRain4Url.toURI());
		testReaderValue.file = testFile.getAbsolutePath();
		testReaderValue.idfield = "ID";
		testReaderValue.tStart = "2022-12-06 17:00";
		testReaderValue.tTimestep = 60;
		testReaderValue.fileNovalue = "-9999";
		testReaderValue.initProcess();

		OmsTimeSeriesIteratorReader predictedFromRReaderValue = new OmsTimeSeriesIteratorReader();
		URL testRainFromR = this.getClass().getClassLoader().getResource("h_fromR_no_trend.csv");
		File testFileFromR = new File(testRainFromR.toURI());
		predictedFromRReaderValue.file = testFileFromR.getAbsolutePath();
		predictedFromRReaderValue.idfield = "ID";
		predictedFromRReaderValue.tStart = "2022-12-06 17:00";
		predictedFromRReaderValue.tTimestep = 60;
		predictedFromRReaderValue.fileNovalue = "-9999";
		predictedFromRReaderValue.initProcess();

		String fId = "id";
		Kriging kriging = new Kriging();
		URL testGridUrl = this.getClass().getClassLoader().getResource("test.shp");
		File testGridFile = new File(testGridUrl.toURI());
		OmsShapefileFeatureReader testReader = new OmsShapefileFeatureReader();
		testReader.file = testGridFile.getAbsolutePath();
		testReader.readFeatureCollection();
		SimpleFeatureCollection testFC = testReader.geodata;

		kriging.inInterpolate = testFC;
		kriging.pSemivariogramType = "exponential";
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
		kriging.fInterpolateid = stationIdField;
		kriging.inHValuesPath = observedFile.getAbsolutePath();
		kriging.cutoffDivide = 20;
//		kriging.nugget = 0;
//		kriging.sill = 15292.38;
//		kriging.range = 82946.36;

		while (reader.doProcess) {
			try {
				reader.nextRecord();
				testReaderValue.nextRecord();
				HashMap<Integer, double[]> id2ValueMap = reader.outData;
				kriging.inData = id2ValueMap;
				kriging.executeKriging();
				predictedFromRReaderValue.nextRecord();
				HashMap<Integer, double[]> predictedGstatR = predictedFromRReaderValue.outData;
				HashMap<Integer, double[]> result = kriging.outData;
				Set<Integer> pointsToInterpolateResult = result.keySet();
				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
				while (iteratorTest.hasNext()) {
					int id = iteratorTest.next();
					double[] values = result.get(id);
					double[] actual = predictedGstatR.get(id);
					assertEquals(actual[0], values[0], 30);
					System.out.println("actual is:" + actual[0] + " evaluate " + values[0]);
				}
			} catch (IOException e) {
				e.printStackTrace();
			} catch (CQLException e) {
				e.printStackTrace();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// writer.close();
	}

	@Test
	public void testKrigingSic97WithTrend() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "id";
		// 100 station to training model
		URL stazioniGridUrl = getClass().getClassLoader().getResource("observed.shp");
		File stazioniGridFile = new File(stazioniGridUrl.toURI());
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		URL observedRain4Url = this.getClass().getClassLoader().getResource("observed_H.csv");
		File observedFile = new File(observedRain4Url.toURI());
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = observedFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2022-12-06 17:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2000-01-01 00:00";
		reader.fileNovalue = "-9999";
		reader.initProcess();
		OmsTimeSeriesIteratorReader testReaderValue = new OmsTimeSeriesIteratorReader();
		URL testRain4Url = this.getClass().getClassLoader().getResource("test_H.csv");
		File testFile = new File(testRain4Url.toURI());
		testReaderValue.file = testFile.getAbsolutePath();
		testReaderValue.idfield = "ID";
		testReaderValue.tStart = "2022-12-06 17:00";
		testReaderValue.tTimestep = 60;
		testReaderValue.fileNovalue = "-9999";
		testReaderValue.initProcess();

		OmsTimeSeriesIteratorReader predictedFromRReaderValue = new OmsTimeSeriesIteratorReader();
		URL testRainFromR = this.getClass().getClassLoader().getResource("h_fromR_no_trend.csv");
		File testFileFromR = new File(testRainFromR.toURI());
		predictedFromRReaderValue.file = testFileFromR.getAbsolutePath();
		predictedFromRReaderValue.idfield = "ID";
		predictedFromRReaderValue.tStart = "2022-12-06 17:00";
		predictedFromRReaderValue.tTimestep = 60;
		predictedFromRReaderValue.fileNovalue = "-9999";
		predictedFromRReaderValue.initProcess();

		String fId = "id";
		Kriging kriging = new Kriging();
		URL testGridUrl = this.getClass().getClassLoader().getResource("test.shp");
		File testGridFile = new File(testGridUrl.toURI());
		OmsShapefileFeatureReader testReader = new OmsShapefileFeatureReader();
		testReader.file = testGridFile.getAbsolutePath();
		testReader.readFeatureCollection();
		SimpleFeatureCollection testFC = testReader.geodata;

		kriging.inInterpolate = testFC;
		kriging.pSemivariogramType = "exponential";
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
		kriging.fInterpolateid = stationIdField;
		kriging.inHValuesPath = observedFile.getAbsolutePath();
		kriging.cutoffDivide = 20;
		kriging.doDetrended = true;
		kriging.fPointZ = "z1";
		kriging.fStationsZ = "z1";

//		kriging.nugget = 0;
//		kriging.sill = 15292.38;
//		kriging.range = 82946.36;

		while (reader.doProcess) {
			try {
				reader.nextRecord();
				testReaderValue.nextRecord();
				HashMap<Integer, double[]> id2ValueMap = reader.outData;
				kriging.inData = id2ValueMap;
				kriging.executeKriging();
				predictedFromRReaderValue.nextRecord();
				HashMap<Integer, double[]> predictedGstatR = predictedFromRReaderValue.outData;
				HashMap<Integer, double[]> result = kriging.outData;
				Set<Integer> pointsToInterpolateResult = result.keySet();
				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
				while (iteratorTest.hasNext()) {
					int id = iteratorTest.next();
					double[] values = result.get(id);
					double[] actual = predictedGstatR.get(id);
					assertEquals(actual[0], values[0], 30);
					System.out.println("actual is:" + actual[0] + " evaluate " + values[0]);
				}
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

}
