Krigings
========

``@author`` Marialaura Bancheri

``@author`` Francesco Serafin (sidereus3), francesco.serafin.3@gmail.com

``@author`` Daniele Andreis, daniele.andreis@gmail.com

``@copyright`` GNU Public License v3 AboutHydrology (Riccardo Rigon)

Description
-----------

The **Krigings** project is designed to perform spatial interpolation
using kriging methods. The implementation includes the following
components:

-  Variogram modeling: Supports multiple types of variograms (e.g.,
   spherical, exponential, Gaussian).
-  Ordinary kriging (OK) and Detrended Kriging for interpolation.
-  Cross-validation module: Evaluates model accuracy through
   leave-one-out cross-validation.
-  Input/Output: Reads input data from shapefile and csv and outputs
   results as CSV/GeoTIFF.

The project is developed using **Java** and built with **Gradle** to
ensure modularity and flexibility.

First version from [@Formetta2014-zn] and [@Bancheri2018-x]

A **built version of the lates release** is available at the `GitHub
release
section <https://github.com/geoframecomponents/Krigings/releases>`__.

Repository and Compatibility
----------------------------

-  The project can be imported into *Eclipse*.
-  It is available on GitHub and can be cloned or downloaded from:
   `GitHub
   Repository <https://github.com/geoframecomponents/Krigings.git>`__
-  **Travis CI** is set up for automated builds and testing, running on
   every commit but only deploying the JAR file when a tagged release is
   pushed to GitHub.
-  The project uses Java ?? (Oracle JDK ??) for compatibility.

Developers’ documentation
-------------------------

This project is a Gradle-based Java project.

.. figure:: doc/ReadMe/gradle.png
   :alt: Gradle logo

   Gradle logo

Building the Project
~~~~~~~~~~~~~~~~~~~~

To build the project, run the following command in the working
directory:

::

   ~ $ gradle build

in the working directory. You will find the built ``.jar`` in
``build/libs``. Once the build process is complete, the generated
``.jar`` file will be located in the ``build/libs`` directory.

Building the ReadMe File
~~~~~~~~~~~~~~~~~~~~~~~~

To convert the ReadMe file from Markdown (``.md``) to RestructuredText
(``.rst``), run:

::

   pandoc doc/ReadMe/ReadMe.md -o ReadMe.rst --bibliography=doc/ReadMe/biblio.bib

This command ensures proper formatting and inclusion of references from
the ``biblio.bib`` file.

Package Structure
~~~~~~~~~~~~~~~~~

**Base Package:** ``org.geoframe.blogpost.kriging``

--------------

``org.geoframe.blogpost.kriging``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Main Class:** - ``Kriging.java`` → Main execution class for kriging
interpolation

--------------

``org.geoframe.blogpost.kriging.interpolationdata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Provides interpolation data, including raster and
vector input handling. **Main Classes:** -
``InterpolationDataProvider.java`` → Interface for providing
interpolation data - ``RasterInterpolationProvider.java`` → Handles
raster-based interpolation - ``VectorInterpolationProvider.java`` →
Handles vector-based interpolation

--------------

``org.geoframe.blogpost.kriging.linearsystemsolver``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Implements solvers for linear systems required in
kriging computations. **Main Class:** -
``SimpleLinearSystemSolverFactory.java`` → Factory class for linear
system solvers

--------------

``org.geoframe.blogpost.kriging.loo`` (Leave-One-Out Kriging)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Description:** Implements Leave-One-Out cross-validation (LOO-CV) for
assessing kriging accuracy. **Main Class:** -
``LeaveOneOutKrigings.java`` → Performs leave-one-out kriging validation

--------------

``org.geoframe.blogpost.kriging.pointcase``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Implements kriging for point-based interpolation.
**Main Class:** - ``KrigingPointCase.java`` → Performs kriging
interpolation on point data, actual implementation of the Kriging
abstract class.

--------------

``org.geoframe.blogpost.kriging.primarylocation``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Manages station selection, residual evaluation, and
spatial processing. **Main Classes:** - ``MaxDistance.java`` → Computes
the maximum interpolation distance - ``Model.java`` → Defines station
selection models - ``NumberOfStations.java`` → Determines the number of
stations used in kriging - ``ResidualsEvaluator.java`` → Compute a
regression and evaluates residuals. - ``StationProcessor.java`` →
Preprocesses station data (filter with StationSelection and compute the
residuals.) - ``StationsSelection.java`` → Implements station selection
strategies

--------------

``org.geoframe.blogpost.kriging.rastercase``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Implements kriging for raster datasets. **Main Class:**
- ``KrigingRasterCase.java`` → Applies kriging on raster data, actual
implementation of Kriging abstract class.

--------------

``org.geoframe.blogpost.kriging.variogram``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Implements variogram modeling, both experimental and
theoretical.

``org.geoframe.blogpost.kriging.variogram.experimental``
                                                        

-  ``ExperimentalVariogram.java`` → Computes the experimental variogram
   from data.

``org.geoframe.blogpost.kriging.variogram.theoretical.curvefitter``
                                                                   

-  ``KrigingParamValidator.java`` → Validates kriging parameters
-  ``VariogramFitter.java`` → Fits variogram models to experimental
   data. **Note**: Initial guess values are set here.
-  ``VariogramFunction.java`` → Implements mathematical functions for
   variogram calculations. **Note**:There is no if statement for
   logarithmic and exponential functions in the case where the x value
   is equal to 0.

``org.geoframe.blogpost.kriging.variogram.theoretical.model``
                                                             

**Supported theoretical variogram models:** - ``Bessel.java`` -
``Circular.java`` - ``Exponential.java`` - ``Gaussian.java`` -
``Spherical.java`` - ``Logarithmic.java`` - ``Pentashperical.java`` -
``Spline.java``

--------------

``org.geoframe.blogpost.kriging.utilities``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Utility functions and helper classes for kriging
computations, including logarithmic transformation and the extraction of
coordinates from a feature collection.

**Note**: This module might be restructured. The logarithmic
transformation could be moved to a separate module to run before and
after (see `issue
#4 <https://github.com/bubbobne/Krigings/issues/4>`__). Additionally,
getCoordinate might be better placed in a different class. —

Linkers’ documentation
----------------------

Integration with OMS3/GEOFrame:

::

   * The project is structured to be compatible with **GEOFrame/OMS3** models.
   * Input data can be retrieved from existing hydrological simulations for further processing.

Users’ documentation
--------------------

(To be completed, it will be published on a Notion static page or on the
blog.)

Under Development
-----------------

PArallelization
~~~~~~~~~~~~~~~

A parallelization method is now under test. If the flag
**parallelComputation** is set to true, a method that uses
**ParallelStream** is employed.

**Note**: This can be useful if you have a large number of points/maps
and are not using tree parallelization.

Future Improvements
-------------------

References
----------

.. container::
   :name: refs
