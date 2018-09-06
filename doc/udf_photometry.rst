.. _udf_photometry:

Photometry of the MUSE UDF fields
=================================

The fitting algorithms in the `imphot` module were tested using MUSE
observations of 10 MUSE fields in the Hubble UDF, denoted UDF01 to
UDF10. The image fitting algorithm worked well for MUSE fields UDF02,
UDF03, UDF08, UDF09, UDF10, where there are no bright point
sources. However bright point sources in fields UDF01, UDF04, UDF05,
UDF06 and UDF07 caused problems until they were excluded. This was
due to two effects.

1. The stars in the UDF fields have sufficient proper motion that they
   had moved significantly in the 11 years that elapsed between the
   original HST UDF observations and the MUSE UDF observations.

2. The HST WCS suffered from a CTE charge loss effect, and this had a
   much larger effect on point sources against a dark background, than
   other sources. The result is that if one tries to scale the fluxes
   of an HST image to match a MUSE image, a different scale factor is
   needed for unresolved point sources, such as stars, than for the
   resolved sources, such as galaxies. Stars need to be scaled by 10%
   to 20% more than galaxies.

Both of these issues break the assumption that the only differences
between an HST and Hubble image of a field are caused by instrumental
differences and atmospheric conditions. To resolve this problem, stars
and QSOs were excluded from the fits, by passing the fitting procedure
a ds9 region file that excludes small circular regions centered on
each star.

The following pages discuss the procedures used to fit the photometric
parameters of each of the MUSE UDF fields. The goal of these pages is
to guide the processing of future observations of these fields, and to
provide real-world examples for processing other fields.

The UDF01 field contains both a star and a QSO, so it provides
examples of the problems noted above.

.. toctree::
   :maxdepth: 2

   udf01.rst
   udf02.rst
   udf03.rst
   udf04.rst
   udf05.rst
   udf06.rst
   udf07.rst
   udf08.rst
   udf09.rst
   udf10.rst
