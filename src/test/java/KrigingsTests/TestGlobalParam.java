package KrigingsTests;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.junit.Test;
import org.opengis.filter.Filter;

import krigingsPointCase.GlobalParameterEvaluator;

public class TestGlobalParam {
	@Test
	public void testGlobalParam() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "ID";
		// 100 station to training model
		URL stationShpUrl = this.getClass().getClassLoader().getResource("MeteoStations.shp");
		File stazioniGridFile = new File("/home/andreisd/Downloads/problemakriging/MeteoStations.shp");
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

	//	URL observedRain4Url = this.getClass().getClassLoader().getResource("precipitation_cleaned12h.csv");
	//	File observedFile = new File("/home/andreisd/Documents/project/GEOFRAME/Krigings/resources/Input/krigings/sic97/precipitation_cleaned12h.csv");
		GlobalParameterEvaluator kriging = new GlobalParameterEvaluator();
		kriging.pSemivariogramType = null;
		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;
	//	kriging.inHValuesPath = observedFile.getAbsolutePath();
		kriging.cutoffDivide = 10;
		kriging.doDetrended = true;
		kriging.fStationsZ = "Altitud";
		kriging.tStart = "1970-10-01 12:00";
		kriging.tTimeStep = 1440;
//		kriging.inNumCloserStations = 5;
		kriging.doLogarithmic = false;
		kriging.inHValuesPath ="/home/andreisd/Downloads/problemakriging/Precipitation.csv";

		kriging.execute();

		// writer.close();
	}
}
