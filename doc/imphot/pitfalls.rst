.. _pitfalls:

Avoiding common pitfalls
========================

Always exclude stars
--------------------

As a general rule, any stars within a field should be masked out when
using the photometric image-fitting method of the `imphot`
module. This can be done using a region file, as described in the
:ref:`regions` section). The reasons for excluding stars are described
below.

The image fitting procedures assume that the differences between a
MUSE image and an HST image, are due to their different PSFs,
calibration errors and pointing errors. The proper motion of a bright
star breaks this assumption if the elapsed time between the MUSE
observation and the HST observation is sufficient for the star to move
more than about a tenth of a pixel.

In the case of the MUSE UDF fields, the HST UDF observations were
performed at the end of 2003, whereas the MUSE UDF observations
started at the end of 2014. In [#f2]_, Pirkal et al (2005), reported
the proper-motions of the following bright unresolved objects in the
MUSE UDF fields (as well as others that are outside the MUSE UDF
fields).

+---------+----------+--------+-------+-----------+------+----------------+
| RA (deg)| Dec (deg)| HST ID | Field | F606W flux| Type | PM (mas/year)  |
+=========+==========+========+=======+===========+======+================+
|53.1630  | -27.7672 | 9397   | UDF01 |  49.0     | QSO  |  0.70 +/- 1.47 |
+---------+----------+--------+-------+-----------+------+----------------+
|53.1580  | -27.7692 | 9230   | UDF01 |  84.0     | Star | 10.26 +/- 0.41 |
+---------+----------+--------+-------+-----------+------+----------------+
|53.1485  | -27.7702 | 9212   | UDF04 |   4.4     | Star | 12.35 +/- 0.54 |
+---------+----------+--------+-------+-----------+------+----------------+
|53.1583  | -27.7949 | 3166   | UDF05 | 193.1     | Star |  8.84 +/- 0.66 |
+---------+----------+--------+-------+-----------+------+----------------+
|53.1766  | -27.7997 | 2150   | UDF06 |  70.5     | Star | 25.16 +/- 0.45 |
+---------+----------+--------+-------+-----------+------+----------------+
|53.1323  | -27.7829 | 5921   | UDF07 | 219.8     | Star |  3.27 +/- 1.26 |
+---------+----------+--------+-------+-----------+------+----------------+

These stars have proper motions from 3.3 to 25.2 mas/year, so during
the 11 years between the HST observations and the MUSE observations,
the slowest of these stars moved 0.04 arcsec, and the fastest moved by
0.3 arcsec. These are significant distances compared to the MUSE pixel
size of 0.2 arcsec, so including these stars in the fitting procedure
will degrade the fit.

Exclude bright QSOs
-------------------

Bright QSOs (Quasi-Stellar-Objects) are point-sources that look
similar to stars. Unlike most stars, they don't have measurable proper
motions, but their fluxes are often significantly variable on both
short and long timescales.  Any image that contains a variable source
will yield a poor fit when the photometric image-fitting method is
applied to it. This is because variable sources break the assumption
that any differences between a MUSE image and a corresponding HST
image, are solely caused by instrumental and atmospheric effects.

In principle the photometric fitting methods can be used if the fitted
region of the image is restricted using a region file. However for
reasons that are not yet understood, fits to bright point sources
yield fitted flux scale factors that are about 10 to 20% higher than
the scale factors of other parts of an image. The other fitted
parameters are more likely to match the parameter fitted to other
parts of the image, but plate distortions or other localized effects
at the position of the QSO make it risky to assume this. In general it
is better to use a region file to mask out the QSO and perform a
global image fit, in order to obtain fitted parameters represent the
average photometric characteristics of the image.

Exclude variable sources
------------------------

For the reasons already discussed in the above sections, any sources
with variable fluxes or positions should be excluded using ds9 region
files, both when aligning images before summing, and when fitting for
the PSF of an image. Weak variable sources may not have much affect on
the fits, but it is easy to mask them out, so this is worth
trying. The quickest way to find out if there are any problematic
sources of this form, is to perform a photometry fit without
specifying a region file, request that the residual image be displayed
(eg. using the :ref:`--display<plotting_options>` option of
:ref:`fit_photometry<fit_photometry>`, or by passing the
``display=True`` to the `~muse_analysis.imphot.fit_image_photometry()`
function), then see if anything significant remains in the residual
image of the fit.

.. [#f2] *Stars in the Hubble Ultra Deep Field*, Pirzkal, N.;
         Sahu, K. C.; Burgasser, A.; Moustakas, L. A.; Xu, C.;
         Malhotra, S.; Rhoads, J. E.; Koekemoer, A. M.; Nelan, E. P.;
         Windhorst, R. A.; Panagia, N.; Gronwall, C.; Pasquali, A.;
         Walsh, J. R., Ap. J. (2005) 622, 319
