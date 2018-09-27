"""
imphot
======

The imphot module provides functions and scripts that can be used to
determine the photometric parameters of individual MUSE images. These
parameters include the width of the point-spread function, the scaling
and zero offset of the flux in the image and the pointing error of the
observation. The photometric parameters of a particular MUSE image are
estimated by comparing that image to an HST image of the same region
of sky. Two algorithms are provided:

1. The `fit_image_photometry()` function uses the HST image as a good
   estimate of the flux distribution on the sky, and sees what
   observing conditions would be needed to reduce the resolution of
   the HST image, and change its scaling and offsets to match the MUSE
   observation.

   To efficiently process many images of a single area, a
   multi-process iterator is also provided that hands the job of
   calling `fit_image_photometry()` on multiple images to multiple
   processes. This is called `FitImagePhotometryMP()`.

2. The `fit_star_photometry()` function fits a Moffat PSF to a given
   star or a bright point source in both the MUSE image and the HST
   image. The PSF parameters are determined directly from the star fit
   in the MUSE image. The flux scaling, flux offset and pointing
   errors are determined by comparing the total fluxes, flux offsets
   and positions of the star fits in the two images.

   To efficiently process many images of a single area, a
   multi-process iterator is also provided that hands the job of
   calling `fit_star_photometry()` on multiple images to multiple
   processes. This is called `FitStarPhotometryMP()`.

In regions that have one or more bright stars, the two methods can be
used to cross check each other. For this purpose, the
`fit_image_photometry()` function can be asked to restrict the area of
the image that it fits to be the same area that is used to fit the
star in `fit_star_photometry()`. This is requested by passing an
optional star=[ra,dec,radius] argument. When both
`fit_image_photometry()` and `fit_star_photometry()` are passed the
same values for this star argument, they operate on exactly the same
pixels of the images. A multiprocessing iterator is also provided for
this case, called `FitImageAndStarPhotometryMP()`.  For each MUSE
image this returns the results from both methods.

Beware that some bright stars suffer from an effect that acts like
saturation. It isn't clear whether this is in the HST observations or
the MUSE observations. The result is that the flux scale-factor that
is found on these stars is higher than elsewhere in the images. Note
that this is seen regardless of whether the fitting is performed by
`fit_image_photometry()` or `fit_star_photometry()`. This implies that
the problem is not an artefact of the fitting method. If
`fit_image_photometry()` is applied to an image that contains a bright
star like this, a poor fit is obtained for both the star and other
objects in the image. The residual image shows that the star ends up
being under-subtracted, while the dimmer objects in the image become
over-subtracted. Again, this indicates that the calibration factor for
the star is different from other sources in the image. In such cases
it may be better to exclude the star from the fit. This can be done by
using the optional ``regions=<filename>`` argument of
`fit_image_photometry()` to specify a ds9 region file that excludes a
small region around the star.

In `fit_image_photometry()`, ds9 region files can be used to mask out
unwanted areas of an image or to selectively indicate which areas of
the image to include in the fit. Note that in ds9, each region can be
marked as inclusive or exclusive, using the Property menu of the
dialog that appears when one double clicks on a region.

The main() functions of the ``fit_photometry``, ``regrid_hst_to_muse`` and
``make_wideband_image`` scripts are also provided. For example, the
``fit_photometry`` script is comprised of the following few lines that
just invoke the main function within the imphot module::

  #!/usr/bin/env python
  import sys
  import imphot
  imphot.fit_photometry_main(sys.argv)

"""

from __future__ import division, absolute_import, print_function

__all__ = ['ds9regions']

# Make all public functions and their documentation visible at the
# package level. In the following, for each file add its public
# functions to __all__, to ensure that their documentation is visible
# when one types help(imphot).  Then delete access to each subpackage,
# to prevent anybody from writing code that depends on the current
# package layout, or on internal functions that may be removed or
# modified in the future.

from .version import __version__

from . import core
from .core import *
__all__.extend(core.__all__)
del(core)

from . import fitimage
from .fitimage import *
__all__.extend(fitimage.__all__)
del(fitimage)

from . import fitstar
from .fitstar import *
__all__.extend(fitstar.__all__)
del(fitstar)

from . import fitboth
from .fitboth import *
__all__.extend(fitboth.__all__)
del(fitboth)

from . import makeimage
from .makeimage import *
__all__.extend(makeimage.__all__)
del(makeimage)

from . import programs
from .programs import *
__all__.extend(programs.__all__)
del(programs)

# Turn the ds9regions class into a sub-module of imphot.

from . import ds9regions

# Remove access to all modules that only contain internal functions.

del mp
