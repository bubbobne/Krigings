package KrigingsTests;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.SchemaException;
import org.geotools.filter.text.cql2.CQLException;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.junit.Test;

import krigingsPointCase.GlobalParameterEvaluator;

public class TestGlobalParam {
	@Test
	public void testGlobalParam() throws URISyntaxException, SchemaException, CQLException, IOException {
		//
		String stationIdField = "ID_shp";
		// 100 station to training model
		String baseUrl = "/home/andreisd/Documents/project/checkKrihingMartin/";
		File stazioniGridFile = new File(baseUrl+"/data/Meteo/Meteo_rain_stations_sel_rem.shp");
		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = stazioniGridFile.getAbsolutePath();
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

	//	URL observedRain4Url = this.getClass().getClassLoader().getResource("precipitation_cleaned12h.csv");
	//	File observedFile = new File("/home/andreisd/Documents/project/GEOFRAME/Krigings/resources/Input/krigings/sic97/precipitation_cleaned12h.csv");
		GlobalParameterEvaluator gbEval = new GlobalParameterEvaluator();
		gbEval.pSemivariogramType = null;
		gbEval.inStations = stationsFC;
		gbEval.fStationsid = stationIdField;
	//	kriging.inHValuesPath = observedFile.getAbsolutePath();
		gbEval.cutoffDivide = 15;
		gbEval.doDetrended = true;
		gbEval.fStationsZ = "elev";
		gbEval.tStart = "2010-01-01 00:00";
		gbEval.tTimeStep = 60;
//		kriging.inNumCloserStations = 5;
		gbEval.doLogarithmic = false;
		gbEval.inHValuesPath = baseUrl+"/data/Meteo//Rain_OMS_format_10_sel_rem.csv";

		gbEval.execute();

		// writer.close();
	}
}
