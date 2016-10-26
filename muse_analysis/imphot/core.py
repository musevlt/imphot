from os.path import basename
from inspect import getargspec
import astropy.units as u
import string
import argparse
import numpy as np
from scipy.ndimage.interpolation import affine_transform
from textwrap import dedent

__all__ = ['FittedValue', 'UserError', 'WorkerError',
           'FIT_IMAGE', 'FIT_STAR', 'FIT_BOTH',
           'ImphotArgumentParser', 'extract_function_args',
           'FittedPhotometry', 'HstFilterInfo', 'rescale_hst_like_muse',
           'regrid_hst_like_muse', 'image_grids_aligned',
           'apply_corrections']

# Set default parameters for the HST Moffat PSF.

_default_hst_fwhm = 0.085
_default_hst_beta = 1.6

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
       True if the value of this parameter was fixed during the\n
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

        if ((FIT_IMAGE | FIT_STAR) & functions):

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

        if FIT_IMAGE & functions:

            self.add_argument('--hst_fwhm', default=_default_hst_fwhm,
                          type=float, metavar="arcsec",
                          help=dedent('''\
                          DEFAULT=%(default)s (arcseconds)
                          The FWHM of a Moffat model of the effective PSF of
                          the HST. The default value that is used if this
                          parameter is not specified, came from Moffat fits to
                          stars within HST UDF images. To obtain the closest
                          estimate to the dithered instrumental PSF, these fits
                          were made to images with the smallest available pixel
                          size (30mas).'''))

            self.add_argument('--hst_beta', default=_default_hst_beta,
                          type=float, metavar="value",
                          help=dedent('''\
                          DEFAULT=%(default)s
                          The beta parameter of a Moffat model of the effective
                          PSF of the HST.  The default value that is used if
                          this parameter is not specified, came from Moffat
                          fits to stars within HST UDF images, as described
                          above for the hst_fwhm parameter. This term is
                          covariant with other terms in the star fits, so there
                          was significant scatter in the fitted values. From
                          this range, a value was selected that yielded the
                          least scatter in the fitted MUSE PSFs in many MUSE
                          images from different MUSE fields and at different
                          wavelengths.'''))

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

            self.add_argument('--taper',  default=9, type=int,
                          metavar="pixels",
                          help=dedent('''\
                          DEFAULT=%(default)s (pixels)
                          This argument controls how transitions
                          between unmasked and masked regions are
                          softened. Because the fitting algorithm
                          replaces masked pixels with zeros, bright
                          sources that are truncted by masked regions
                          cause sharp changes in brightness that look
                          like real features and bias the fitted
                          position error. To reduce this effect,
                          pixels close to the boundary of a masked
                          region are tapered towards zero over a
                          distance specified by the --taper argument.
                          The method used to smooth the transition
                          requires that the --taper argument be an odd
                          number of pixels, so if an even-valued
                          integer is specified, this is quietly
                          rounded up to the next highest odd
                          number. Alternatively, the softening
                          algorithm can be disabled by specifying 0
                          (or any value below 2).'''))

            self.add_argument('--extramask',  default=None, type=_str_or_none,
                          metavar="text",
                          help=dedent('''\
                          DEFAULT=%(default)s
                          If the value of this argument is not the
                          word "none", then it should name a FITS file
                          that contains a mask image to be combined
                          with the mask of the MUSE image.
                          Specifically, this FITS file should have an
                          IMAGE extension called 'DATA' and the image
                          in that extension should have the same
                          dimensions and WCS coordinates as the MUSE
                          images. The pixels of the image should be
                          integers, with 0 used to denote unmasked
                          pixels, and 1 used to denote masked pixels.'''))


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

            self.add_argument('--apply', action='store_true',
                          help=dedent('''\
                          Derive corrections from the fitted position errors
                          and calibration errors, apply these to the MUSE
                          image, and write the resulting image to a FITS
                          file. The name of the output file is based on the
                          name of the input file, by replacing its ".fits"
                          extension with "_aligned.fits".  If the input muse
                          image was not read from a file, then a file called
                          "muse_aligned.fits" is written in the current
                          directory. Also see the --resample option.'''))

            self.add_argument('--resample', action='store_true',
                          help=dedent('''\
                          By default the --apply option corrects position
                          errors by changing the coordinate reference pixel
                          (CRPIX1,CRPIX2) without changing any pixel values.
                          Alternatively, this option can be used to shift
                          the image by resampling its pixels, without
                          changing the coordinates of the pixels.'''))

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
    muse : `mpdaf.obj.Image`
       The MUSE image that the fit was performed on.
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
    rms_error : float
       The root-mean square of the pixel residuals of the fit in the
       MUSE image, in the same units as the pixels of the original
       MUSE image.

    Attributes
    ----------
    name : str
       The basename of the MUSE FITS without the .fits extension.
    fit_report : str
       A printable report on the fit from the least-squares\n
       fitting function.
    scale : `FittedValue`
       The best-fit value and error of the calibration scale\n
       factor, MUSE/HST.
    bg : `FittedValue`
       The best-fit value and error of the calibration offset,\n
       MUSE-HST.
    dx : `FittedValue`
       The best-fit value and error of the x-axis pointing offset,\n
       MUSE.x-HST.x (arcsec). This is the distance that features\n
       in the HST image had to be moved to the right, in the\n
       direction of increasing x-axis pixel index, to line them up\n
       with the same features in the MUSE image. One way to\n
       correct the pointing error of the MUSE observation, is to\n
       divide 'dx' by the pixel size along the x-axis (usually\n
       0.2 arcsec), then add the resulting pixel offset to the\n
       CRPIX1 header parameter.
    dy : `FittedValue`
       The best-fit value and error of the y-axis pointing offset,\n
       MUSE.y-HST.y (arcsec). This is the distance that features\n
       in the HST image had to be moved upwards, in the direction\n
       of increasing y-axis pixel index, to line them up with the\n
       same features in the MUSE image. One way to correct the\n
       pointing error of the MUSE observation, is to divide 'dy'\n
       by the pixel size along the y-axis (usually 0.2 arcsec),\n
       then add the resulting pixel offset to the CRPIX2 header\n
       parameter.
    dra : `FittedValue`
       The right-ascension error (arcsec) that corresponds to the\n
       pointing error dx,dy. This is the angular distance that\n
       features in the HST image had to be moved towards increased\n
       right-ascension, to line them up with the same feaures in\n
       the MUSE image. One way to correct the pointing error of\n
       the MUSE observation is to subtract 'dra' from the CRVAL1\n
       header value of the MUSE observation.
    ddec : `FittedValue`
       The declination error (arcsec) that corresponds to the\n
       pointing error dx,dy. This is the angular distance that\n
       features in the HST image had to be moved towards increased\n
       declination, to line them up with the same feaures in\n
       the MUSE image. One way to correct the pointing error of\n
       the MUSE observation is to subtract 'ddec' from the CRVAL2\n
       header value of the MUSE observation.
    fwhm : `FittedValue`
       The best-fit value and error of the FWHM of the Moffat\n
       PSF.
    beta : `FittedValue`
       The best-fit value and error of the beta parameter of the\n
       Moffat PSF.
    rms_error : float
       The root-mean square of the pixel residuals of the fit in\n
       the MUSE image, in the same units as the pixels of the\n
       original MUSE image.

    """

    def __init__(self, method, muse, fit_report, scale, bg, dx, dy,
                 fwhm, beta, rms_error):
        self.method = method
        self.name = basename(muse.filename).replace(".fits","")
        self.fit_report = fit_report
        self.scale = scale
        self.bg = bg
        self.dx = dx
        self.dy = dy
        self.fwhm = fwhm
        self.beta = beta
        self.rms_error = rms_error

        # Convert the X and Y axis corrections, and their standard
        # deviations, from arcsec to pixel counts.

        steps = muse.wcs.get_step(unit=u.arcsec)
        p = np.array([dy.value, dx.value]) / steps
        s = np.array([dy.stdev, dx.stdev]) / steps

        # We need to Right Ascension and declination offsets that
        # correspond, approximately, to the dx,dy pointing errors.
        # Calculate these by adding the offsets to the pixel at the
        # center of the image.

        center_pix = np.array(np.asarray(muse.shape)/2.0)
        center_dec_ra = muse.wcs.pix2sky(center_pix, unit=u.arcsec)[0]
        ddec,dra = muse.wcs.pix2sky(center_pix + p, unit=u.arcsec)[0] -\
                   center_dec_ra
        sdec,sra = muse.wcs.pix2sky(center_pix + s, unit=u.arcsec)[0] -\
                   center_dec_ra

        # Record the right ascension and declination corrections in
        # arcseconds.

        self.dra = FittedValue(value = dra, stdev = sra,
                               fixed = dx.fixed and dy.fixed)
        self.ddec = FittedValue(value = ddec, stdev = sdec,
                                fixed = dx.fixed and dy.fixed)

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
            name_width = len(self.name)
            s = "#%-*s Method  Flux    FWHM    beta     Flux    x-offset  y-offset  ra-offset dec-offset  RMS\n" % (name_width-1, " MUSE observation ID")
            s += "#%*s         scale  arcsec           offset    arcsec    arcsec    arcsec    arcsec    error\n" % (name_width-1, "")
            s +=  "#%s ------ ------- ------- ------- --------- --------- --------- --------- --------- -------\n" % ('-' * (name_width - 1))

        # Format the fitted values.

        s += "%s %6s % 7.4f % 7.4f % 7.4f % 9.5f % 9.5f % 9.5f % 9.5f % 9.5f %7.4f" % (
            self.name, self.method, self.scale.value, self.fwhm.value,
            self.beta.value, self.bg.value,
            self.dx.value, self.dy.value, self.dra.value, self.ddec.value,
            self.rms_error)
        return s

class HstFilterInfo(object):
    """An object that contains the filter characteristics of an HST
    image.

    Parameters
    ----------
    hst : `mpdaf.obj.Image` or str
       The HST image to be characterized, or the name of an HST filter,
       such as "F606W.

    Attributes
    ----------
    filter_name : str
       The name of the filter.
    abmag_zero : float
       The AB magnitude that corresponds to the zero electrons/s\n
       from the camera (see ``https://archive.stsci.edu/prepds/xdf``).
    photflam : float
       The calibration factor to convert pixel values in\n
       electrons/s to erg cm-2 s-1 Angstrom-1.\n
       Calculated using::

         photflam = 10**(-abmag_zero/2.5 - 2*log10(photplam) - 0.9632)

       which is a rearranged form of the following equation from\n
       ``http://www.stsci.edu/hst/acs/analysis/zeropoints``::

         abmag_zero = -2.5 log10(photflam)-5 log10(photplam)-2.408

    photplam : float
       The pivot wavelength of the filter (Angstrom)\n
       (from ``http://www.stsci.edu/hst/acs/analysis/bandwidths``)
    photbw : float
       The effective bandwidth of the filter (Angstrom)\n
       (from ``http://www.stsci.edu/hst/acs/analysis/bandwidths``)

    """

    # Create a class-level dictionary of HST filter characteristics.
    # See the documentation of the attributes for the source of these
    # numbers.

    _filters = {
        "F606W" :  {"abmag_zero" : 26.51,    "photplam"  : 5921.1,
                    "photflam"   : 7.73e-20, "photbw"    : 672.3},
        "F775W" :  {"abmag_zero" : 25.69,    "photplam"  : 7692.4,
                    "photflam"   : 9.74e-20, "photbw"    : 434.4},
        "F814W" :  {"abmag_zero" : 25.94,    "photplam"  : 8057.0,
                    "photflam"   : 7.05e-20, "photbw"    : 652.0},
        "F850LP" : {"abmag_zero" : 24.87,    "photplam"  : 9033.1,
                    "photflam"   : 1.50e-19, "photbw"    : 525.7}
    }

    def __init__(self, hst):

        # If an image has been given, get the name of the HST filter
        # from the FITS header of the HST image, and convert the name
        # to upper case. Otherwise get it from the specified filter
        # name.

        if isinstance(hst, str):
            self.filter_name = hst.upper()
        elif "FILTER" in hst.primary_header:
            self.filter_name = hst.primary_header['FILTER'].upper()
        elif ("FILTER1" in hst.primary_header and
              hst.primary_header['FILTER1'] != 'CLEAR1L'):
            self.filter_name = hst.primary_header['FILTER1'].upper()
        elif ("FILTER2" in hst.primary_header and
              hst.primary_header['FILTER2'] != 'CLEAR2L'):
            self.filter_name = hst.primary_header['FILTER2'].upper()
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

def rescale_hst_like_muse(hst, muse, inplace=True):
    """Rescale an HST image to have the same flux units as a given MUSE image.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       The HST image to be resampled.
    muse : `mpdaf.obj.Image` or `mpdaf.obj.Cube`
       A MUSE image or cube with the target flux units.
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
    """Resample an HST image onto the spatial coordinate grid of a given
    MUSE image or MUSE cube.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       The HST image to be resampled.
    muse : `mpdaf.obj.Image` of `mpdaf.obj.Cube`
       The MUSE image or cube to use as the template for the HST image.
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

    # If a MUSE cube was provided, extract a single-plane image to use
    # as the template.

    if muse.ndim > 2:
        muse = muse[0,:,:]

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

def apply_corrections(im, imfit, corrections="xy,scale,zero", resample=False,
                      inplace=False):

    """Correct an image for pointing-errors, calibration errors and/or
    zero offsets, using the fitted errors returned by a preceding call
    to `imphot.fit_image_photometry()` or `imphot.fit_star_photometry()`.

    Note that this function is called on the user's behalf by
    `imphot.fit_image_photometry()` and `imphot.fit_star_photometry()` if
    those functions are given the argument, apply=True.

    Parameters
    ----------
    im : `mpdaf.obj.Image`
       The image to be corrected.
    imfit : `imphot.FittedPhotometry`
       Fitted image properties returned by
       `imphot.fit_image_photometry()` or
       `imphot.fit_star_photometry()`.
    corrections : str
       A comma-separated list of the corrections that are to be
       applied, indicated by the following words.

         xy     - The x-axis and y-axis pointing error.
         scale - The flux scale-factor.
         zero  - The flux zero-offset.

       The default string is "x,y,scale,zero", which applies
       all of the corrections.
    resample : bool
       This parameter controls how pointing errors are corrected.
       They can be corrected either by changing the coordinates of the
       pixels (resample=False), or by resampling the pixels to shift
       the image within the original coordinate grid (resample=True).
    inplace : bool
       By default a new image is returned, and the input image is
       not changed. However by setting inplace=True, the image is
       corrected within the input container.

    Returns
    -------
    out : `mpdaf.obj.Image`
       The corrected image.

    """

    # Get a container for the output image.

    out = im if inplace else im.copy()

    # Split the list of corrections to be applied.

    if corrections is None:
        pars = []
    else:
        pars = corrections.split(",")

    # Correct the pointing errors?

    if "xy" in pars:

        # Get the size of the image pixels in arcseconds.

        pixh, pixw = im.get_step(unit=u.arcsec)

        # Convert the corrections from arcseconds to pixels.

        dx = imfit.dx.value / pixw
        dy = imfit.dy.value / pixh

        # Resample the image to correct its pointing errors?

        if resample:

            # Since we only need to shift the image, the affine transform
            # matrix is the identity matrix, and only the pixel offsets need
            # to be changed.

            affine_matrix = np.identity(2)
            affine_offset = np.array([dy, dx])

            # Resample the MUSE image to correct its pointing error.

            data = affine_transform(out.data.filled(0.0), affine_matrix,
                                    affine_offset, prefilter=False, order=2)
            if out.var is None:
                var = None
            else:
                var = affine_transform(out.var.filled(0.0), affine_matrix,
                                       affine_offset, prefilter=False, order=2)
            if out.mask is np.ma.nomask:
                mask = np.ma.nomask
            else:
                mask = affine_transform(out.mask.astype(float), affine_matrix,
                                  affine_offset, prefilter=False, order=2) > 0.1

            # Install the modified arrays.

            out._data = data
            out._var = var
            out._mask = mask

        # Change the coordinates of the image to correct for position errors?

        else:

            out.wcs.set_crpix1(out.wcs.get_crpix1() + dx)
            out.wcs.set_crpix2(out.wcs.get_crpix2() + dy)


    # Correct the zero offset of the fluxes?

    if "zero" in pars:
        out._data -= imfit.bg.value

    # Correct the image for flux scaling errors?

    if "scale" in pars:
        out._data /= imfit.scale.value
        out._var /= imfit.scale.value**2

    # Return the corrected image.

    return out
