package KrigingsTests;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.junit.Test;
import org.opengis.filter.Filter;
import org.opengis.filter.FilterFactory2;
import org.opengis.filter.expression.Expression;

import oldPointCase.Krigings;

public class TestOldKriging {

	@Test
	public void testKrigingSic97WithTrend() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "ID";
		// 100 station to training model
		URL stationShpUrl = this.getClass().getClassLoader().getResource("MeteoStations.shp");
		File stazioniGridFile = new File(stationShpUrl.toURI());
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		URL observedRain4Url = this.getClass().getClassLoader().getResource("precipitation_normal_score.csv");
		File observedFile = new File(observedRain4Url.toURI());
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = observedFile.getAbsolutePath();
		reader.idfield = "ID";
		reader.tStart = "2017-12-12 08:00";
		reader.tTimestep = 60;
		// reader.tEnd = "2017-12-12 05:00";
		reader.fileNovalue = "-9999.0";
		reader.initProcess();

		String fId = "ID";
		Krigings kriging = new Krigings();

		int idToCheck = 25;
		Filter filter = CQL.toFilter("ID = " + idToCheck);

		kriging.inInterpolate = stationsFC.subCollection(filter);
		kriging.pSemivariogramType = "exponential";
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
		kriging.fInterpolateid = stationIdField;
		kriging.doDetrended = true;
		kriging.fPointZ = "z_dem";
		kriging.fStationsZ = "z_dem";
//			kriging.inNumCloserStations = 5;
		kriging.nugget = 0.0;
		kriging.sill = 0.3;
		kriging.range = 6000;
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = "resources/Output/krigings/PointCase/stations_normal_" + idToCheck + ".csv";
		writer.tStart = reader.tStart;
		writer.tTimestep = reader.tTimestep;
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
					// assertEquals(actual[0], values[0], 2);
					System.out.println("actual is:" + actual[0] + " evaluate " + values[0]);
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
}
