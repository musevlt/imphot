from os.path import basename
from inspect import getargspec
import astropy.units as u
import string
import argparse
import numpy as np
from textwrap import dedent

__all__ = ['FittedValue', 'UserError', 'WorkerError',
           'FIT_IMAGE', 'FIT_STAR', 'FIT_BOTH',
           'ImphotArgumentParser', 'extract_function_args',
           'FittedPhotometry', 'HstFilterInfo', 'rescale_hst_like_muse',
           'regrid_hst_like_muse', 'image_grids_aligned']

# Define a class that will hold the fitted value of a single parameter.

class FittedValue(object):
    """A parameter value returned by a least-squares fitting function.

    Parameters
    ----------
    param : `lmfit.parameter.Parameter`
        The best-fit value of a parameter.

    Attributes
    ----------
    value : float
       The best-fit value.
    stdev : float
       The estimated 1-sigma uncertainty in the best-fit value.
    fixed : bool
       True if the value of this parameter was fixed during the
       fit. False if the parameter was free to vary.
    """

    def __init__(self, param=None, value=None, stdev=None, fixed=None):
        if param is not None:
            self.value = param.value
            self.stdev = param.stderr
            self.fixed = not param.vary
        else:
            self.value = value
            self.stdev = stdev
            self.fixed = fixed


# Define an exception for errors found in the input files or other data.

class UserError(Exception):
    """An exception that is raised when an error occurs that is due to
    user-specified inputs, rather than a programming error."""

    def __init__(self, message):
        self.message = message
    def __str__(self):
        return self.message

class WorkerError(Exception):
    """An exception class for reporting exceptions that occur
    in separate worker processes.

    Parameters
    ----------
    message : str
       A short message about the exception.
    traceback : list of str
       The traceback of the exception, in the form of a list
       of text lines.

    """

    def __init__(self, message, traceback):
        self.message = message
        self.traceback = string.join(traceback)

    def __str__(self):
        return self.message + "\n" + self.traceback

def _string_to_float_or_none(s):
    """Parse a command-line argument string that can contain either
    a floating point value, or the word "none".

    If the string contains a valid floating point value, that value is
    returned. If instead it contains the word "none" or "None", or "",
    then None is returned. If it contains neither a floating point
    value nor the word none, then an `UserError` exception is
    raised.

    Parameters
    ----------
    s  : str
       The string that is expected to contain a floating point
       value or the word "none".
    Returns
    -------
    out : float or None
       The value of a floating point number, or None.

    """
    try:
        return None if s.lower() in ("none", "") else float(s)
    except ValueError as e:
        raise UserError(e.message)

def _str_or_none(s):
    """Parse a command-line argument string that can contain either
    a string or the word "none".

    If the string contains the word none, return None.

    Parameters
    ----------
    s  : str
       The string to be parsed.
    Returns
    -------
    out : str or None
       The parsed string.

    """
    return None if s.lower() == "none" else s

# Enumerate the two photometry fitting functions, and their combination,
# using separate bits for FIT_IMAGE and FIT_STAR.

FIT_IMAGE = 1                     # bit 0
FIT_STAR = 2                      # bit 1
FIT_BOTH = FIT_IMAGE + FIT_STAR   # bits 0 and 1

# Describe the sub-set of command-line arguments that are
# useful for calling fit_image_photometry().

class ImphotArgumentParser(argparse.ArgumentParser):
    """A parser of command-line arguments, that recognizes options that
    are used by `fit_image_photometry()` and/or
    `fit_star_photometry()`.

    For example, to create a parser that recognizes the optional arguments
    of both `fit_image_photometry()` and `fit_star_photometry()`,
    one would do the following::

      parser = imphot.ImphotArgumentParser(imphot.FIT_IMAGE +
                                     imphot.FIT_STAR, prog=argv[0])

    To use this parser to process the command-line arguments of a program,
    one would then use the parse_args() function of argparse.ArgumentParser
    as follows::

      options = parser.parse_args(argv[1:])

    This returns an object that contains initialized variables for
    each of the arguments recognized by `imphot.ImphotArgumentParser`. From
    this one can obtain separate dictionaries of the keyword/value
    pairs of the arguments of the `fit_image_photometry()` and
    `fit_star_photometry()` functions by using the
    `extract_function_args()` function as follows::

      imfit_kwargs = imphot.extract_function_args(options,
                                 imphot.fit_image_photometry)
      stfit_kwargs = imphot.extract_function_args(options,
                                 imphot.fit_star_photometry))

    These dictionaries can then be passed as keyword/value arguments
    to `fit_image_photometry()` and `fit_star_photometry()` as follows::

      imfit = imphot.fit_image_photometry(hst, muse, **imfit_kwargs)
      stfit = imphot.fit_star_photometry(hst, muse, **stfit_kwargs)

    Parameters
    ----------
    functions : int
       The fitting functions that we want to obtain arguments
       for, expressed as a sum of one or both of
       `imphot.FIT_IMAGE` and `imphot.FIT_STAR`.

    """

    def __init__(self, functions, *args, **kwargs):

        # Instantiate the super-class OptionParser object.

        super(ImphotArgumentParser, self).__init__(*args, formatter_class=argparse.RawTextHelpFormatter, **kwargs)

        # Define command-line arguments according to the function
        # that is being configured.

        if FIT_IMAGE & functions:
            self.add_argument('--regions',  default="star", type=_str_or_none,
                          metavar="file-star-notstar-or-none",
                          help=dedent('''\
                          DEFAULT=%(default)s

                          This can be "none", to indicate that all pixels
                          should be used in the global image fit, "star" to
                          restrict the fit to pixels within any optional circle
                          defined by the --star argument, "notstar" to restrict
                          the image fit to pixels outside any circle defined by
                          the --star argument, or the filename of a ds9 region
                          file.

                          This option can be used to exclude damaged regions of
                          an image, or to exclude sources that would degrade
                          the global PSF fit, such as saturated stars, stars
                          with significant proper motion, or variable sources.

                          Alternatively this option can also be used to
                          restrict the fit to one or more objects, by masking
                          everything except small regions around these objects.

                          Only ds9 circle, ellipse and box regions are
                          supported. Other types of ds9 regions and
                          configuration parameters, are simply ignored. For
                          each region, be careful to tell ds9 whether you want
                          the region to be included or excluded. Also be
                          careful to specify either fk5 or physical pixel
                          coordinates.

                          If the value of this argument is "star", and the
                          --star argument has also been provided, then a region
                          is substituted that only includes the circular region
                          specified by the --star argument.

                          Similarly, if "notstar" is passed, and the --star
                          argument has been provided, then a region is
                          substituted that excludes the circular region set by
                          the --star argument.'''))

        if ((FIT_IMAGE | FIT_STAR) & functions):

            self.add_argument('--star', nargs=3, default=None, type=float,
                          metavar=("ra", "dec", "radius"),
                          help=dedent('''\
                          DEFAULT=%(default)s
                          Perform photometry fits to a star at a specified
                          position. This is done in addition to the global
                          image fitting procedure, so two sets of
                          photometry results are reported when this option
                          is selected.

                          The option expected 3 values. These are the right
                          ascension and declination of the star in decimal
                          degrees, and the radius of the area over which to
                          perform the fit, in arcseconds.

                          If the --regions argument is also present, and
                          its value is the word "star", then the imaging
                          fitting procedure is also restricted to the
                          region of this star. This can be used as a
                          consistency check, as both methods should yield
                          similar results.'''))

        if FIT_IMAGE & functions:

            self.add_argument('--segment', action='store_true',
                          help=dedent('''\
                          Ignore areas that don't contain significant objects
                          by ignoring pixels that are below the median value in
                          a morphologically opened version of the HST image.'''))

            self.add_argument('--fix_scale',  default=None,
                          metavar="factor",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s
                          Use this option to fix the calibration scale
                          factor, (MUSE_flux / HST_flux) to the specified
                          value while fitting. The default value is "none",
                          which means that the parameter will be fitted.'''))

            self.add_argument('--fix_bg',  default=None,
                          metavar="offset",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s
                          Use this option to fix the calibration zero-offset
                          (MUSE_flux - HST_flux) to the specified value while
                          fitting. The default value is "none", which means
                          that the parameter will be fitted.'''))

            self.add_argument('--fix_dx',  default=None,
                          metavar="arcsec",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s (arcseconds)
                          Use this option to fix the x-axis pointing offset,
                          (MUSE_x - HST_x) to the specified value while
                          fitting. The default value is "none", which means
                          that the parameter will be fitted.'''))

            self.add_argument('--fix_dy',  default=None,
                          metavar="arcsec",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s (arcseconds)
                          Use this option to fix the y-axis pointing offset,
                          (MUSE_y - HST_y) to the specified value while
                          fitting. The default value is "none", which means
                          that the parameter will be fitted.'''))

            self.add_argument('--fix_fwhm',  default=None,
                          metavar="arcsec",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s (arcseconds)
                          Use this option to fix the FWHM of the Moffat PSF
                          to the specified value while fitting. The default
                          value is "none", which means that the parameter
                          will be fitted.'''))

            self.add_argument('--fix_beta',  default=2.5,
                          metavar="value",
                          type=_string_to_float_or_none,
                          help=dedent('''\
                          DEFAULT=%(default)s
                          Use this option to fix the beta exponent of the
                          Moffat PSF to the specified value while fitting.
                          The default value is 2.5. Change this to "none"
                          if you wish this parameter to be fitted.'''))

            self.add_argument('--margin',  default=2.0, type=float,
                          metavar="arcsec",
                          help=dedent('''\
                          DEFAULT=%(default)s (arcseconds)
                          The width of a margin of zeros to add around the
                          image before processing. A margin is needed because
                          most of the processing is performed using discrete
                          Fourier transforms, which are periodic in the width
                          of the image. Without a margin, features at one edge
                          of the image would spill over to the opposite edge of
                          the image when a position shift was applied, or when
                          features were widened by convolving them with a
                          larger PSF. The margin width should be the maximum of
                          the largest expected position error between the two
                          input images, and the largest expected PSF width.'''))

            self.add_argument('--save', action='store_true',
                          help=dedent('''\
                          Save the result images of each input image to FITS
                          files.'''))

        # Describe command-line arguments that are common to
        # both fit_image_photometry() and fit_star_photometry().


        if ((FIT_IMAGE | FIT_STAR) & functions):

            self.add_argument('--display', action='store_true',
                          help=dedent('''\
                          Display the images, FFTs and star fits, if any.'''))

            self.add_argument('--nowait', action='store_true',
                          help=dedent('''\
                          Don't wait for the user to interact with each
                          displayed plot before continuing.'''))

            self.add_argument('--hardcopy', default="none",
                          type=_str_or_none, metavar="format-or-none",
                          help=dedent('''\
                          Write hardcopy plots of the fitting results to files
                          that have the specified graphics format (eg. "pdf",
                          "jpg", "png", "eps").  Plots of the fits will be
                          written to filenames that start with the name of the
                          MUSE input file, after removing any .fits suffix,
                          followed by either "_image_fit.<suffix>" for the plot
                          of the image fit, or "_star_fit.<suffix>" for plots
                          of any star fits.'''))

            self.add_argument('--title',  default=None, type=_str_or_none,
                          metavar="text",
                          help=dedent('''\
                          DEFAULT=%(default)s
                          Either a plot title, "none" to request the default
                          title, or "" to request that no title be displayed
                          above the plots.'''))


def extract_function_args(options, function):
    """Given a dictionary of key/value pairs or a Namespace object
    returned by `argparse.ArgumentParser.parse_args()`, extract the
    subset of argument values that are recognized by a specified
    function and return them in a dictionary of keyword/value
    pairs. The caller can then pass this dictionary to the
    specified function as keyword/value arguments.

    Parameters
    ----------
    options : `argparse.Namespace` or `dict`
       The return value of a call to `argparse.ArgumentParser.parse_args()`.
    function : function
       The function whose arguments are needed.

    Returns
    -------
    out : dict
       A dictionary of keyword/value pairs that can be passed
       to the specified function.
    """
    # Convert the options object into a dictionary if necessary.
    if isinstance(options, dict):
        opts = options
    else:
        opts = vars(options)

    # Get the argument list of the function.

    argspec = getargspec(function)

    # Get the intersection of the two dictionaries.

    return {key: opts[key] for key in argspec.args if key in opts}

# Define a class that holds fitted photometry parameters.

class FittedPhotometry(object):
    """The superclass of `FittedImagePhotometry` and `FittedStarPhotometry`
    which contains the fitted photometry parameters of a MUSE image,
    the identity of the original MUSE image, and a report from the
    least-squares fit.

    Parameters
    ----------
    method : str
        A short word that identifies the method used to fit the photometry.
    filename : str
        The name of the original MUSE file.
    fit_report : str
        A multi-line string that contains a verbose report about the
        final least-squares fit or fits.
    scale : `FittedValue`
        The best-fit value and error of the calibration scale
        factor, MUSE/HST.
    bg : `FittedValue`
        The best-fit value and error of the calibration offset,
        MUSE-HST.
    dx : `FittedValue`
        The best-fit value and error of the x-axis pointing offset,
        MUSE.x-HST.x.
    dy : `FittedValue`
        The best-fit value and error of the y-axis pointing offset,
        MUSE.y-HST.y.
    fwhm : `FittedValue`
        The best-fit value and error of the FWHM of the Moffat PSF.
    beta : `FittedValue`
        The best-fit value and error of the beta parameter of the
        Moffat PSF.

    Attributes
    ----------
    name : str
        The basename of the MUSE FITS without the .fits extension.
    fit_report : str
        A printable report on the fit from the least-squares
        fitting function.
    scale : `FittedValue`
        The best-fit value and error of the calibration scale
        factor, MUSE/HST.
    bg : `FittedValue`
        The best-fit value and error of the calibration offset,
        MUSE-HST.
    dx : `FittedValue`
        The best-fit value and error of the x-axis pointing offset,
        MUSE.x-HST.x.
    dy : `FittedValue`
        The best-fit value and error of the y-axis pointing offset,
        MUSE.y-HST.y.
    fwhm : `FittedValue`
        The best-fit value and error of the FWHM of the Moffat PSF.
    beta : `FittedValue`
        The best-fit value and error of the beta parameter of the
        Moffat PSF.

    """

    def __init__(self, method, filename, fit_report, scale, bg, dx, dy,
                 fwhm, beta):
        self.method = method
        self.name = basename(filename).replace(".fits","")
        self.fit_report = fit_report
        self.scale = scale
        self.bg = bg
        self.dx = dx
        self.dy = dy
        self.fwhm = fwhm
        self.beta = beta

    def __str__(self):
        return "Report of HST %s photometry fit of MUSE observation %s\n" % (self.method, self.name) + self.fit_report

    def summary(self, header=True):
        """Return a string that summarizes the results, optionally
        listed under a 3-line heading of column descriptions.

        Parameters
        ----------
        heading : bool
           True to include 3 header lines above the values (default=True)
        """

        s = ""

        # Include header lines?

        if header:
            name_width = 34
            s = "# MUSE observation ID              Method    Flux    FWHM    beta      Flux  x-offset  y-offset\n"
            s += "#                                           scale     (\")            offset       (\")       (\")\n"
            s += "#--------------------------------- ------  ------  ------  ------  --------  --------  --------\n"

        # Format the fitted values.

        s += "%34s %6s % 6.4f % 7.4f % 7.4f % 9.5f % 9.5f % 9.5f" % (
            self.name, self.method, self.scale.value, self.fwhm.value,
            self.beta.value, self.bg.value,
            self.dx.value, self.dy.value)
        return s

class HstFilterInfo(object):
    """An object that contains the filter characteristics of an HST
    image.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       The HST image to be characterized.

    Attributes
    ----------
    filter_name : str
       The name of the filter.
    abmag_zero : float
       The AB magnitude that corresponds to the zero electrons/s
       from the camera. This comes from (from ``https://archive.stsci.edu/prepds/xdf``)
    photflam : float
       The calibration factor to convert pixel values in
       electrons/s to erg cm-2 s-1 Angstrom-1.
       Calculated using::

         photflam = 10**(-abmag_zero/2.5 - 2*log10(photplam) - 0.9632)

       which is a rearranged form of the following equation from
       ``http://www.stsci.edu/hst/acs/analysis/zeropoints``::

         abmag_zero = -2.5 log10(photflam)-5 log10(photplam)-2.408

    photplam : float
       The pivot wavelength of the filter (Angstrom)
       (from ``http://www.stsci.edu/hst/acs/analysis/bandwidths``)
    photbw : float
       The effective bandwidth of the filter (Angstrom)
       (from ``http://www.stsci.edu/hst/acs/analysis/bandwidths``)
    psf_fwhm : float
      The full-width at half maximum of a Moffat profile fitted to
      a star in the HST XDF (arcsec)
    psf_beta : float
      The beta value of a Moffat profile fitted to a star in the
      HST XDF.

    """

    # Create a class-level dictionary of HST filter characteristics.
    # See the documentation of the attributes for the source of these
    # numbers.

    _filters = {
        "F606W" :  {"abmag_zero" : 26.51,    "photplam"  : 5921.1,
                    "photflam"   : 7.73e-20, "photbw"    : 672.3,
                    "psf_fwhm"   : 0.085,    "psf_beta"  : 1.6},
        "F775W" :  {"abmag_zero" : 25.69,    "photplam"  : 7692.4,
                    "photflam"   : 9.74e-20, "photbw"    : 434.4,
                    "psf_fwhm"   : 0.085,    "psf_beta"  : 1.6},
        "F814W" :  {"abmag_zero" : 25.94,    "photplam"  : 8057.0,
                    "photflam"   : 7.05e-20, "photbw"    : 652.0,
                    "psf_fwhm"   : 0.085,    "psf_beta"  : 1.6},
        "F850LP" : {"abmag_zero" : 24.87,    "photplam"  : 9033.1,
                    "photflam"   : 1.50e-19, "photbw"    : 525.7,
                    "psf_fwhm"   : 0.085,    "psf_beta"  : 1.6}
    }

    def __init__(self, hst):

        # Get the name of the HST filter from the FITS header of the
        # HST image, and convert the name to upper case.

        if "FILTER" in hst.primary_header:
            self.filter_name = hst.primary_header['FILTER'].upper()
        else:
            raise UserError("Missing FILTER keyword in HST image header")

        # Get the dictionary of the characteristics of the filter.

        if self.filter_name in self.__class__._filters:
            info = self.__class__._filters[self.filter_name]
        else:
            raise UserError("Unrecognized filter (%s) in HST image header" %
                            self.filter_name)

        # Record the characteristics of the filer.

        self.abmag_zero = info['abmag_zero']
        self.photflam = info['photflam']
        self.photplam = info['photplam']
        self.photbw = info['photbw']
        self.psf_fwhm = info['psf_fwhm']
        self.psf_beta = info['psf_beta']

        for attr in info:
            setattr(self, attr, info[attr])

def rescale_hst_like_muse(hst, muse, inplace=True):
    """Rescale an HST image to have the same flux units as a given MUSE image.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       The HST image to be resampled.
    muse : `mpdaf.obj.Image`
       The MUSE image with the target flux units.
    inplace : bool
       (This defaults to True, because HST images tend to be large)
       If True, replace the contents of the input HST image object
       with the rescaled image.
       If False, return a new Image object that contains the rescaled
       image.

    Returns
    -------
    out : `mpdaf.obj.Image`
       The rescaled HST image.

    """

    # Operate on a copy of the input image?

    if not inplace:
        hst = hst.copy()

    # Get the characteristics of the HST filter.

    filt = HstFilterInfo(hst)

    # Calculate the calibration factor needed to convert from
    # electrons/s in the HST image to MUSE flux-density units.

    cal = filt.photflam * u.Unit("erg cm-2 s-1 Angstrom-1").to(muse.unit)

    # Rescale the HST image to have the same units as the MUSE image.

    hst.data *= cal
    if hst.var is not None:
        hst.var *= cal**2
    hst.unit = muse.unit

    return hst

def regrid_hst_like_muse(hst, muse, inplace=True):
    """Resample an HST image onto the coordinate grid of a given MUSE image.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       The HST image to be resampled.
    muse : `mpdaf.obj.Image`
       The MUSE image to use as the template for the HST image.
    inplace : bool
       (This defaults to True, because HST images tend to be large)
       If True, replace the contents of the input HST image object
       with the resampled image.
       If False, return a new Image object that contains the resampled
       image.

    Returns
    -------
    out : `mpdaf.obj.Image`
       The resampled HST image.

    """

    # Operate on a copy of the input image?

    if not inplace:
        hst = hst.copy()

    # Mask the zero-valued blank margins of the HST image.

    np.ma.masked_inside(hst.data, -1e-10, 1e-10, copy=False)

    # Resample the HST image onto the same coordinate grid as the MUSE
    # image.

    return hst.align_with_image(muse, cutoff=0.0, flux=True, inplace=True)

def image_grids_aligned(im1, im2, tolerance=0.01):
    """Return True if two images sample the same array of sky positions
    to within a specified tolerance.

    Parameters
    ----------
    im1  :  mpdaf.obj.Image
       The first of the images to be compared.
    im2  :  mpdaf.obj.Image
       The second of the images to be compared.
    tolerance : float
       The tolerance within which positions must match, in arcseconds.

    Returns
    -------
    out : bool
       True if the images have matching coordinate grids; False if they
       don't.

    """

    # Convert the tolerance from arcseconds to degrees.

    atol = tolerance / 3600.0

    # Get the dimensions of the first image.

    ny,nx = im1.shape

    # Get the index of the central pixel.

    center = (ny//2,nx//2)

    # Are the coordinate conversion matrices of the two images the same?

    if not np.allclose(im1.wcs.get_cd(), im2.wcs.get_cd(),
                       rtol=0.0, atol=atol):
        return False

    # Are the array dimensions of the images the same?

    if not np.allclose(im1.shape, im2.shape):
        return False

    # Are the coordinates of the central pixel of each image the same?

    if not np.allclose(im1.wcs.pix2sky(center), im2.wcs.pix2sky(center),
                       rtol=0.0, atol=atol):
        return False

    # The images have matching coordinate grids.

    return True
