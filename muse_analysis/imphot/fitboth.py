from .mp import _FitPhotometryMP
from .fitimage import fit_image_photometry
from .fitstar import fit_star_photometry
from .core import (extract_function_args, UserError, _default_hst_fwhm,
                   _default_hst_beta)
from mpdaf.obj import Image

__all__ = ['fit_image_and_star_photometry', 'FitImageAndStarPhotometryMP']

def fit_image_and_star_photometry(hst, muse, star, regions="star",
                                  fix_scale=None, fix_bg=None, fix_dx=None,
                                  fix_dy=None, fix_fwhm=None, fix_beta=None,
                                  hst_fwhm=_default_hst_fwhm,
                                  hst_beta=_default_hst_beta,
                                  margin=2.0, segment=False, display=False,
                                  nowait=False, hardcopy=None, title=None,
                                  save=False, fig=None):

    """Given a MUSE image and an HST image that has been regridded and
    aligned onto the same coordinate grid as the MUSE image, use the
    HST image as a calibrated reference to fit for the flux-scale, the
    FWHM and beta parameter of a Moffat PSF, the position offset, and
    the background offset of the MUSE image relative to the HST
    image. The fit is performed two ways; first by calling
    `fit_image_photometry()`, then again by calling
    `fit_star_photometry()`. Both sets of results are returned.

    Parameters
    ----------
    hst : `mpdaf.obj.Image` or filename
       An HST image of the same area as the MUSE image, and that has been
       gridded onto the same pixel grid as the MUSE image. This can be
       given as an MPDAF Image object, or by the filename of the FITS file.
    muse : `mpdaf.obj.Image` or filename
       The MUSE image to be characterized. This can be
       given as an MPDAF Image object, or by the filename of the FITS file.
    star : (float, float, float)
       This argument indicates the position of a single star within the
       MUSE image and the radius of a circular area of pixels around it
       that should be used in the fitting process by `fit_star_photometry()`.

       If the optional regions argument is set to the special string,
       "star", then the star's position and the associated radius are
       also passed to `fit_image_photometry()` to indicate that the image
       fit should also be restricted to the same circular area around the
       star.

       The value of the star argument is a tuple of the right
       ascension and declination of the star in decimal degrees, and
       the radius of the area over which to perform the fit, in
       arcseconds.
    regions : str or iterable or None
       This can be None, to indicate that no regions are needed in the
       global image fit, "star" to restrict the image fit to pixels
       within any region defined by the star=(ra,dec,radius) argument,
       "notstar" to restrict the image fit to pixels outside any region
       defined the by star=(ra,dec,radius) argument, the name of a
       filename of the ds9 region file, or an iterable that returns
       successive lines of a ds9 region file.

       These regions are passed to `fit_image_photometry()`. They can
       be used to exclude problematic areas of an image, or to exclude
       sources that would degrade the global PSF fit, such as
       saturated stars, stars with significant proper motion, and
       variable sources.  Alternatively it can be used to restrict the
       fit to one or more objects by masking everything except small
       regions around these objects. In particular, if the value of
       the regions argument is "star", then the fit is restricted to
       the circular area specified by star argument.

       Only ds9 circle, ellipse and box regions are supported. Other region
       types and most of the configuration parameters found in ds9 region
       files, are simply ignored. For each region, be careful to tell ds9
       whether you want the region to be included or excluded. Also be
       careful to request either fk5 or physical pixel coordinates, because
       other coordinate systems are not supported.
    fix_scale : float or None
       The calibration scale factor, (MUSE_flux / HST_flux) is fixed
       to the specified value while fitting, unless the value is None.
       This parameter affects the image fit but not the star-profile fit.
    fix_bg : float or None
       The calibration zero-offset, (MUSE_flux - HST_flux) is fixed
       to the specified value while fitting, unless the value is None.
       This parameter affects the image fit but not the star-profile fit.
    fix_dx : float or None
       The x-axis pointing offset, (MUSE_x - HST_x) is fixed
       to the specified value while fitting, unless the value is None.
       This parameter affects the image fit but not the star-profile fit.
    fix_dy : float or None
       The y-axis pointing offset, (MUSE_y - HST_y) is fixed
       to the specified value while fitting, unless the value is None.
       This parameter affects the image fit but not the star-profile fit.
    fix_fwhm : float or None
       The FWHM of the Moffat PSF is fixed to the specified value
       while fitting, unless the value is None.
       This parameter affects both the image fit and the FWHM of
       the star profile that is fit to the MUSE image.
    fix_beta : float or None
       The beta exponent of the Moffat PSF is fixed to the specified value
       while fitting, unless the value is None.
       This parameter affects both the image fit and the beta parameter of
       the star profile that is fit to the MUSE image.
    hst_fwhm : float
       The FWHM of a Moffat model of the effective PSF of the HST.
       The default value that is used if this parameter is not
       specified, came from Moffat fits to stars within HST UDF
       images. To obtain the closest estimate to the dithered
       instrumental PSF, these fits were made to images with the
       smallest available pixel size (30mas).
    hst_beta : float
       The beta parameter of a Moffat model of the effective PSF of
       the HST.  The default value that is used if this parameter is
       not specified, came from Moffat fits to stars within HST UDF
       images, as described above for the hst_fwhm parameter.
    margin : float
       The width (arcsec) of a margin of zeros to add around the image
       before processing. A margin is needed because most of the
       processing is performed using discrete Fourier transforms, which
       are periodic in the width of the image. Without a margin, features
       at one edge of the image would spill over to the opposite edge of
       the image when a position shift was applied, or when features were
       widened by convolving them with a larger PSF. The margin width
       should be the maximum of the largest expected position error
       between the two input images, and the largest expected PSF width.
    segment : bool
       If True, ignore areas that don't contain significant objects by
       ignoring pixels that are below the median value in a
       morphologically opened version of the HST image.
    display : bool
       If True (the default), display the plot.
    hardcopy : str or None
       If this is a non-empty string, then it should contain a
       graphics file suffix supported by matplotlib, such as "pdf",
       "jpg", "png" or "eps". Plots of the fits will be written to
       filenames that start with the name of the MUSE input file,
       after removing any .fits suffix, followed by
       "_image_fit.<suffix>" for the plot of the image fit, and
       "_star_fit.<suffix>" for the plots of the star fits.
    nowait : bool
       When this argument is False, wait for the user to dismiss
       the plot before returning. This allows the user to interact
       with the plot. When this argument is True, the plot will
       dissapear as soon as another plot is drawn and the user will
       not be able to interact with it.
    title : str or None
       A specific plot title, or None to request the default title.
       Specify "" if no title is wanted.
    save : bool
       If True, save the result images of each input image to
       FITS files.

    Returns
    -------
    out : (`FittedImagePhotometry`, `FittedStarPhotometry`)
       An object that contains the fitted parameters from both
       `fit_image_photometry()` and `fit_star_photometry()`.

    """

    # If needed read the FITS files into MPDAF Image objects.

    try:
        if isinstance(hst, str):
            hst = Image(hst)
    except Exception as e:
        raise UserError("Error reading HST file (%s)" % e.message)

    try:
        if isinstance(muse, str):
            muse = Image(muse)
    except Exception as e:
        raise UserError("Error reading MUSE file (%s)" % e.message)

    # Check the star argument.

    try:
        ra, dec, radius = star
    except Exception:
        raise UserError("The star argument should be a tuple of (ra,dec,radius)")
    if radius <= 0.0:
        raise UserError("The radius given to the star argument must be > 0")

    # Perform the image fit.

    image_results = fit_image_photometry(hst, muse, regions=regions,
                                         fix_scale=fix_scale, fix_bg=fix_bg,
                                         fix_dx=fix_dx, fix_dy=fix_dy,
                                         fix_fwhm=fix_fwhm, fix_beta=fix_beta,
                                         hst_fwhm=hst_fwhm, hst_beta=hst_beta,
                                         margin=margin, segment=segment,
                                         display=display,
                                         nowait=nowait,
                                         hardcopy=hardcopy, title=title,
                                         star=star, save=save, fig=fig)

    # Perform the star fit.

    star_results = fit_star_photometry(hst, muse, star, fix_fwhm=fix_fwhm,
                                       fix_beta=fix_beta, display=display,
                                       nowait=nowait, hardcopy=hardcopy,
                                       title=title, fig=fig)

    return (image_results, star_results)

class FitImageAndStarPhotometryMP(_FitPhotometryMP):
    """A multiprocessing iterator that creates a pool or worker processes
    to repeatedly call `fit_image_and_star_photometry()` for each of a
    list of MUSE files, returning successive results from
    `fit_image_and_star_photometry()` via the iterator, as they become
    available.

    Parameters
    ----------
    hst_filename : str
       The name of a FITS fie that contains an HST has the same
       coordinate grid as the MUSE image FITS files that are to be
       processed.
    muse_filenames : list of str
       A list of filenames of the FITS files of the MUSE images that
       are to be processed. These must all be of the same field,
       with the same image coordinate grid as the HST file.
    kwargs : dict
       A dictionary of keyword/value arguments to be passed to
       fit_image_and_star_photometry().
    nworker : int
       The number of worker processes to use. The default is
       0, which creates multiprocessing.cpu_count() processes.
       Alternatively, if a negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.

    """
    def __init__(self, hst_filename, muse_filenames, kwargs, nworker=0):

        # Ensure that kwargs contains the mandatory star argument,
        # and that it contains the expected fields.

        if "star" not in kwargs:
            raise UserError("Missing argument: star=(ra,dec,radius)")
        try:
            ra, dec, radius = kwargs["star"]
        except Exception:
            raise UserError("The star argument should be a tuple of (ra,dec,radius)")
        if radius <= 0.0:
            raise UserError("The radius given to the star argument must be > 0")
        # Instantiate the parent _FitPhotometryMP object.

        super(FitImageAndStarPhotometryMP, self).__init__(
            hst_filename, muse_filenames,
            cmd_fn=fit_image_and_star_photometry, cmd_kwargs=kwargs,
            nworker=nworker)

    def next(self):
        """Return the results from the next image in the list of input MUSE
        images.

        Returns
        -------
        out : `FittedStarPhotometry`
           The fitting results from the next image.

        """

        return super(FitImageAndStarPhotometryMP, self).next()
