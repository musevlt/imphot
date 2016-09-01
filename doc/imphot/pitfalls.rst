.. _pitfalls:

Avoiding common pitfalls
========================

Stars should be excluded
------------------------

As a general rule, any stars within a field should be masked out when
using the photometric image-fitting method of the `imphot`
module. This can be done using a region file, as described in the
:ref:`regions` section). The reasons for excluding stars are
described below.

The image fitting procedure assumes that all differences between a
MUSE image and an HST image, are due to their different PSFs, their
different calibration errors and different pointing errors. Stars can
break this assumption in two ways:

1. Stars that have significant proper motions may appear to have moved
   between the HST and MUSE observations of a field. For example there
   were 11 years between the HST observations of the UDF and the start
   of the MUSE observations of the UDF. In that time, all of the stars
   in the UDF moved significantly. Estimated proper motions of these
   stars were reported in [#f1]_ by Pirkal et al (2005), and the ones
   seen in the MUSE UDF fields are listed below.

   +----------+-----------+--------+-------+--------------------+------+----------------+
   | RA (deg) | Dec (deg) | HST ID | Field | Flux [erg/Å/s/cm²] | Type | PM (mas/year)  |
   +==========+===========+========+=======+====================+======+================+
   | 53.1630  | -27.7672  | 9397   | UDF01 |  12.3e-18          | QSO  |  0.70 ± 1.47   |
   +----------+-----------+--------+-------+--------------------+------+----------------+
   | 53.1580  | -27.7692  | 9230   | UDF01 |  21.0e-18          | Star | 10.26 ± 0.41   |
   +----------+-----------+--------+-------+--------------------+------+----------------+
   | 53.1485  | -27.7702  | 9212   | UDF04 |   1.1e-18          | Star | 12.35 ± 0.54   |
   +----------+-----------+--------+-------+--------------------+------+----------------+
   | 53.1583  | -27.7949  | 3166   | UDF05 |  49.0e-18          | Star |  8.84 ± 0.66   |
   +----------+-----------+--------+-------+--------------------+------+----------------+
   | 53.1766  | -27.7997  | 2150   | UDF06 |  17.3e-18          | Star | 25.16 ± 0.45   |
   +----------+-----------+--------+-------+--------------------+------+----------------+
   | 53.1323  | -27.7829  | 5921   | UDF07 |  56.3e-18          | Star |  3.27 ± 1.26   |
   +----------+-----------+--------+-------+--------------------+------+----------------+

   The proper motions of these stars range from 3.3 to 25.2
   mas/year. Over the 11 years that elapsed between the HST
   observations and the MUSE observations, this corresponded to an
   accumulated motion of 0.04 arcsec to 0.3 arcsec. These are
   significant distances compared to the MUSE pixel size of 0.2
   arcsec, so these stars significantly degrade the image fitting
   procedure if they aren't excluded.

2. Another problem is that the HST WFC suffered from a problem called
   CTE charge loss. Although this affected all sources in an image, it
   is was much worse for faint unresolved point sources against a dark
   background. In [#f2]_, Reiss (2003) determined that for faint stars
   against a background, as much of 10% of the flux of the source
   could be lost. This was before the deeper UDF observations, and for
   the UDF, with its darker backgrounds and sensitivity to fainter
   stars, the effect is worse. The photometric image-fitting technique
   reveals that the flux of stars in the HST UDF images typically need
   to be scaled by 10% to 20% more than other sources in the same
   images, to match the MUSE image. If stars are left in the fit, then
   this results in a poor fit, unless the fit is restricted to a small
   area centered on a single star.

QSOs should also be excluded
----------------------------

Bright QSOs (Quasi-Stellar-Objects) are point-sources that look
similar to stars. Unlike most stars, they don't have measurable proper
motions, but their fluxes are often significantly variable on both
short and long timescales. This, along with the CTE effect which was
described above for stars, results in poor fits when any bright QSOs
are included during the photometric image-fitting procedure.

If the photometric image fitting procedure is limited to a small
region centered on a QSO, then this should yield a good fit, which can
be used to estimate the pointing errors in the image. However the flux
scale that it reports should not be used to scale the rest of the
image. Beware that an estimate of the pointing error that is obtained
from a single QSO, will be affected by any localized geometric
distortions at the position of the QSO in the image. For this reason,
it may be better to use a region file to mask out the QSO and perform
a global image fit on the other sources in the image. This should
yield an estimate of the pointing error that averages out localized
distortions.

Exclude variable sources
------------------------

For the reasons already discussed in the above sections, any sources
with varying fluxes or positions should be excluded using ds9 region
files, both when aligning images before summing, and when fitting for
the PSF of an image. Weak variable sources may not have much affect on
the fits, but it is easy to mask them out, so this is worth
trying. The quickest way to find out if there are any problematic
sources of this form, is to perform a photometry image-fit without
specifying a region file, and request that the residual image be
displayed, using the :ref:`--display<plotting_options>` option of the
:ref:`fit_photometry<fit_photometry>` script. Any problematic sources
will then be visible in the residual image of the fit.

.. [#f1] *Stars in the Hubble Ultra Deep Field*, Pirzkal, N.;
         Sahu, K. C.; Burgasser, A.; Moustakas, L. A.; Xu, C.;
         Malhotra, S.; Rhoads, J. E.; Koekemoer, A. M.; Nelan, E. P.;
         Windhorst, R. A.; Panagia, N.; Gronwall, C.; Pasquali, A.;
         Walsh, J. R., Ap. J. (2005) 622, 319

.. [#f2] *On-orbit Calibration of ACS CTE Corrections for Photometry*,
         Riess, A, *Instrument Science Report ACS 2003-009* (2003).
