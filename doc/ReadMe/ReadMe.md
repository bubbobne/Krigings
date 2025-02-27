# Krigings

`@author` Marialaura Bancheri

`@author` Francesco Serafin (sidereus3), francesco.serafin.3@gmail.com

`@author` Daniele Andreis, daniele.andreis@gmail.com


`@copyright` GNU Public License v3 AboutHydrology (Riccardo Rigon)

## Description

The Krigings project is designed to perform spatial interpolation using kriging methods.
The implementation includes the following components:

  *  Variogram modeling: Supports multiple types of variograms (e.g., spherical, exponential, Gaussian).
  *  Ordinary kriging (OK) and Detrended Kriging for interpolation.
  *  Cross-validation module: Evaluates model accuracy through leave-one-out cross-validation.
  *  Input/Output: Reads input data from shapefile and csv and outputs results as CSV/GeoTIFF.

The project is developed using Java and built with Gradle to ensure modularity and flexibility.

First version from [@Formetta2014-zn] and [@Bancheri2018-x]

A built version of the lates release is available at the [GitHub release section](https://github.com/geoframecomponents/Krigings/releases).

## Documentation

#### Repository and Compatibility

* The project can be imported into Eclipse.
* It is available on GitHub and can be cloned or downloaded from: [GitHub Repository](https://github.com/geoframecomponents/Krigings.git)
* Travis CI is set up for automated builds and testing, running on every commit but only deploying the JAR file when a tagged release is pushed to GitHub.
* The project uses Java ?? (Oracle JDK ??) for compatibility.


### Developers' documentation
This project is a **Gradle-based Java project**.

![Gradle logo](doc/ReadMe/gradle.png)

#### Building the Project

To build the project, run the following command in the working directory:

    ~ $ gradle build

in the working directory. You will find the built `.jar` in `build/libs`.
Once the build process is complete, the generated `.jar` file will be located in the `build/libs` directory.

#### Building the ReadMe File

To convert the ReadMe file from Markdown (`.md`) to RestructuredText (`.rst`), run:

    pandoc doc/ReadMe/ReadMe.md -o ReadMe.rst --bibliography=doc/ReadMe/biblio.bib

This command ensures proper formatting and inclusion of references from the `biblio.bib` file.



#### The Code

**Base Package:** `org.geoframe.blogpost.kriging`

---

## Package Structure

### `interpolationdata`
**Description:** Provides interpolation data, including raster and vector input handling.
**Main Classes:**
- `InterpolationDataProvider.java` → Interface for providing interpolation data
- `RasterInterpolationProvider.java` → Handles raster-based interpolation
- `VectorInterpolationProvider.java` → Handles vector-based interpolation

---

### `linearsystemsolver`
**Description:** Implements solvers for linear systems required in kriging computations.
**Main Class:**
- `SimpleLinearSystemSolverFactory.java` → Factory class for linear system solvers

---

### `loo` (Leave-One-Out Kriging)
**Description:** Implements Leave-One-Out cross-validation (LOO-CV) for assessing kriging accuracy.
**Main Class:**
- `LeaveOneOutKrigings.java` → Performs leave-one-out kriging validation

---

### `pointcase`
**Description:** Implements kriging for point-based interpolation.
**Main Class:**
- `KrigingPointCase.java` → Performs kriging interpolation on point data

---

### `primarylocation`
**Description:** Manages station selection, residual evaluation, and spatial processing.
**Main Classes:**
- `MaxDistance.java` → Computes the maximum interpolation distance
- `Model.java` → Defines station selection models
- `NumberOfStations.java` → Determines the number of stations used in kriging
- `ResidualsEvaluator.java` → Evaluates residuals from kriging results
- `StationProcessor.java` → Preprocesses station data
- `StationsSelection.java` → Implements station selection strategies

---

### `rastercase`
**Description:** Implements kriging for raster datasets.
**Main Class:**
- `KrigingRasterCase.java` → Applies kriging on raster data

---

### `variogram`
**Description:** Implements variogram modeling, both experimental and theoretical.

#### Subpackages

##### `experimental`
- `ExperimentalVariogram.java` → Computes the experimental variogram from data

##### `theoretical.curvefitter`
- `KrigingParamValidator.java` → Validates kriging parameters
- `VariogramFitter.java` → Fits variogram models to experimental data
- `VariogramFunction.java` → Implements mathematical functions for variogram calculations

##### `theoretical.model`
**Supported theoretical variogram models:**
- `Bessel.java`
- `Circular.java`
- `Exponential.java`
- `Gaussian.java`
- `Spherical.java`
- `Logarithmic.java`
- `Pentashperical.java`
- `Spline.java`

---

### `utilities`
**Description:** Utility functions and helper classes for kriging computations.

---

### `org.geoframe.blogpost.kriging`
**Main Class:**
- `Kriging.java` → Main execution class for kriging interpolation

---

### Linkers' documentation

Integration with OMS3/GEOFrame:

    * The project is structured to be compatible with GEOFrame/OMS3 models.
    * Input data can be retrieved from existing hydrological simulations for further processing.


### Users' documentation

(To be completed, it will be published on a Notion static page or on the blog.)


###  Future Improvements



## References

::: {#refs}
:::
