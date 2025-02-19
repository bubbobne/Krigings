package KrigingsTests;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.junit.Test;
import org.opengis.filter.Filter;

public class TestMM {
	@Test
	public void testKrigingMM() throws URISyntaxException, SchemaException, CQLException, IOException {
//		//
//		String stationIdField = "basinid";
//		// 100 station to training model
//		String baseUrl = "/home/andreisd/Documents/project/bacino_po/";
//		File stazioniGridFile = new File(baseUrl+"/data/four_regions_1033.shp");
//		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
//		stationsReader.file = stazioniGridFile.getAbsolutePath();
//		stationsReader.readFeatureCollection();
//		SimpleFeatureCollection stationsFC = stationsReader.geodata;
//
//		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
//		OmsTimeSeriesIteratorWriter w = new OmsTimeSeriesIteratorWriter();
//
//		File observedFile = new File(baseUrl+"/data/four_regions_precipitation.csv");
//		File wFile = new File(baseUrl+"/data/out.csv");
//		reader.file = observedFile.getAbsolutePath();
//		reader.idfield = "ID";
//		reader.tStart = "1991-08-24 00:00";
//		reader.tEnd = "2001-08-25 00:00";
//
//		reader.tTimestep = 1440;
//		// reader.tEnd = "2017-12-12 05:00";
//		reader.fileNovalue = "-9999.0";
//		reader.initProcess();
//		w.file = wFile.getAbsolutePath();
//		w.tStart = "1991-08-24 00:00";
//		w.tTimestep = 1440;
//		// reader.tEnd = "2017-12-12 05:00";
//		w.fileNovalue = "-9999.0";
//		w.inTablename = "ID";
//	   
//		Kriging kriging = new Kriging();
//		File pointGridFile = new File(baseUrl+"2170/centroid_ID_2170.shp");
//		OmsShapefileFeatureReader pointReader = new OmsShapefileFeatureReader();
//		pointReader.file = pointGridFile.getAbsolutePath();
//		pointReader.readFeatureCollection();
//		SimpleFeatureCollection pointFC = pointReader.geodata;
//		kriging.doIncludeZero=true;
//		kriging.inInterpolate = pointFC;
//		//kriging.pSemivariogramType = "exponential";
//		kriging.inStations = stationsFC;
//		kriging.fStationsid = "ID";
//		kriging.fInterpolateid = "basinid";
//        kriging.inHValuesPath= observedFile.getAbsolutePath();
//		kriging.cutoffDivide = 15;
//		kriging.doDetrended = true;
//		kriging.fPointZ = "elev_m";
//		kriging.fStationsZ = "quota";
//		kriging.inNumCloserStations = 15;
//		kriging.doLogarithmic = false;
//		kriging.boundedToZero = true;
//		kriging.doIncludeZero=true;
//        kriging.pAGlobal= 243520;
//        kriging.pSGlobal=122;
//        kriging.pNugGlobal=13;
//        
//       kriging.pAGlobalDeTrended=114000;
//       kriging.pSGlobalDeTrended=66;
//       kriging.pNugGlobalDeTrended=5;
//
//
//        kriging.tTimeStep=1440;
//        kriging.globalDetrendedVariogramType="exponential";
//        kriging.globalVariogramType="spherical";
//		
//		
//		while (reader.doProcess) {
//			try {
//				reader.nextRecord();
//				System.out.println("new T:" + reader.tCurrent);
//				HashMap<Integer, double[]> id2ValueMap = reader.outData;
//
//				HashMap<Integer, double[]> predictedGstatR = new HashMap<Integer, double[]>();
//				kriging.inData = id2ValueMap;
//				kriging.executeKriging();
//				HashMap<Integer, double[]> result = kriging.outData;
//				Set<Integer> pointsToInterpolateResult = result.keySet();
//				Iterator<Integer> iteratorTest = pointsToInterpolateResult.iterator();
//				w.inData = kriging.outData;
//				w.writeNextLine();
//				while (iteratorTest.hasNext()) {
//					int id = iteratorTest.next();
//					double[] values = result.get(id);
//					double[] actual = predictedGstatR.get(id);
//					// assertEquals(actual[0], values[0], 2);
//					System.out.println(" evaluate " + values[0]);
//					HashMap<Integer, double[]> params = kriging.outVariogramParams;
//					System.out.println("n"+params.get(0)[0]);
//					System.out.println("s"+params.get(1)[0]);
//					System.out.println("r"+params.get(2)[0]);
//					
//					System.out.println("d"+params.get(3)[0]);
//					System.out.println("f"+params.get(4)[0]);
//					
//
//				}
//
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
//
//
//		reader.close();
//        w.close();
//		// writer.close();
	}

	
}
