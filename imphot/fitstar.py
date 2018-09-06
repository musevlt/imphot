from os.path import basename
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from lmfit import Model

from mpdaf.obj import Image

from .core import (UserError, FittedValue, FittedPhotometry,
                   image_grids_aligned, apply_corrections)
from .mp import _FitPhotometryMP

__all__ = ['fit_star_photometry', 'FittedStarPhotometry', 'FitStarPhotometryMP']


def fit_star_photometry(hst, muse, star, fix_fwhm=None, fix_beta=None,
                        display=False, nowait=False, hardcopy=None,
                        title=None, fig=None, apply=False, resample=False):
    """Given a MUSE image and an HST image that have been regridded
    and aligned onto the same coordinate grid as the MUSE image, fit
    moffat profiles to a specified star within both of these images,
    then compare the fits to determine the X,Y pointing offset, pixel
    zero-offset and pixel scaling error of the MUSE image. Also use
    the fitted FWHM and beta parameters to indicate the characteristics
    of the MUSE PSF.

    Optionally, if apply=True is passed to this function, pointing
    corrections and calibration corrections are derived from the
    fitted parameters, and a corrected version of the MUSE image is
    written to a new FITS file. See also the resample option, which
    controls how pointing errors are corrected.

    Parameters
    ----------
    hst : `mpdaf.obj.Image`
       An HST image of the same area as the MUSE image, and that has been
       gridded onto the same pixel grid as the MUSE image.
    muse : `mpdaf.obj.Image`
       The MUSE image to be characterized.
    star : float, float, float
       The central right-ascension, declination and radius of the area
       within which to fit the stellar brightness profile. The right
       ascension and declination are in degrees. The radius is in
       arcseconds.
    fix_fwhm : float or None
       The FWHM of the Moffat PSF in the MUSE image is fixed to the
       specified value while fitting, unless the value is None. This
       does not affect the PSF that is fitted to the star in the HST
       image.
    fix_beta : float or None
       The beta exponent of the Moffat PSF in the MUSE image is fixed to
       the specified value while fitting, unless the value is None. This
       does not affect the PSF that is fitted to the star in the HST
       image.
    display : bool
       If True (the default), display the image.
    hardcopy : str or None
       Unless this is an empty string or the word, "none", then it
       should contain a graphics file suffix supported by matplotlib,
       such as "pdf", "jpg", "png" or "eps". Plots of the star fits
       will be written to a filename that starts with the name of the
       MUSE input file, after removing any .fits suffix, followed by
       "_star_fit.<suffix>".
    nowait : bool
       When this argument is False, wait for the user to dismiss
       the plot before returning. This allows the user to interact
       with the plot. When this argument is True, the plot will
       dissapear as soon as another plot is drawn and the user will
       not be able to interact with it.
    title : str or None
       A specific plot title, None or "none" to request the default title,
       or "" if no title is wanted.
    apply : bool
       If True, derive corrections from the fitted position errors and
       calibration errors, apply these to the MUSE image, and write
       the resulting image to a FITS file. The name of the output file
       is based on the name of the input file, by replacing its
       ".fits" extension with "_aligned.fits".  If the input muse
       image was not read from a file, then a file called
       "muse_aligned.fits" is written in the current directory.
    resample : bool
       When apply==True, this argument determines how position errors
       are corrected. If resample=False, then the coordinate reference
       pixel (CRPIX1, CRPIX2) is adjusted to change the coordinates of
       the pixels without changing any pixel values. If resample=True,
       then the pixel values are resampled to shift the image without
       changing the coordinates of the pixels.

    Returns
    -------
    out : `FittedStarPhotometry`
       The results are returned in an object.

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

    # Require that the image arrays of the HST and MUSE images
    # sample the same coordinate grid.

    if not image_grids_aligned(hst, muse):
        raise UserError("The pixels of the HST and MUSE images are not aligned.")

    # Get the central coordinates and radius of the fitting area.

    try:
        ra, dec, radius = star
    except Exception:
        raise UserError("The star argument should be a tuple of (ra,dec,radius)")
    if radius <= 0.0:
        raise UserError("The radius given to the star argument must be > 0")

    # Get the Y-axis and X-axis array dimensions and the corresponding
    # pixel widths.

    ny, nx = muse.shape
    dy, dx = muse.wcs.get_step(unit=u.arcsec)

    # Calculate the area of a pixel.

    pixel_area = dx * dy

    # Get the floating point pixel index of the specified star position.

    ystar, xstar = muse.wcs.sky2pix([dec, ra])[0]

    # Check that the star is within the image.

    if ystar < 0 or ystar >= nx or xstar < 0.0 or xstar >= ny:
        raise UserError("The specified star position is not within the image.")

    # Get the X and Y indexes of each pixel in the image, relative
    # to the specified star position.

    xx, yy = np.meshgrid((np.arange(nx, dtype=float) - xstar) * dx,
                         (np.arange(ny, dtype=float) - ystar) * dy)

    # Get a 2D array of x**2+y**2 for each pixel of the image,
    # relative to the star position.

    rsq = xx**2 + yy**2

    # Find the flattened array indexes of pixels with radii less
    # than the maximum radius that we want to fit.

    indexes = np.nonzero(rsq.ravel() < radius**2)[0]

    # Sort the indexes into order of increasing radius.

    sorted_indexes = indexes[np.argsort(rsq.ravel()[indexes])]

    # Extract the x, y coordinates of the above pixels.

    x = xx.ravel()[sorted_indexes]
    y = yy.ravel()[sorted_indexes]

    # Get the values of the MUSE pixels at each element of x and y.

    muse_values = muse.data.filled(0.0).ravel()[sorted_indexes]

    # Create the 2D Moffat model that will be fit to the star.  Set
    # the initial guesses for dx and dy to the flux-weighted means of
    # x and y, and the initial guess at the moffat flux to the sum of
    # the pixel values.

    fitmod = Model(_circular_moffat_profile,
                   independent_vars=["x", "y", "pixel_area"])
    if fix_fwhm is None:
        fitmod.set_param_hint('fwhm', value=0.5, min=0.0)
    else:
        fitmod.set_param_hint('fwhm', value=fix_fwhm,
                              min=0.0, vary=False)
    if fix_beta is None:
        fitmod.set_param_hint('beta', value=3.0)
    else:
        fitmod.set_param_hint('beta', value=fix_beta, vary=False)
    fitmod.set_param_hint('dx', value=np.average(x, weights=muse_values))
    fitmod.set_param_hint('dy', value=np.average(y, weights=muse_values))
    fitmod.set_param_hint('flux', value=muse_values.sum(), min=0.0)
    fitmod.set_param_hint('bg', value=0.0)
    fitmod.make_params()

    # Fit the model to the MUSE data.

    muse_results = fitmod.fit(muse_values, x=x, y=y, pixel_area=pixel_area,
                              weights=10.0)

    # Calculated the fitted MUSE star model at x,y.

    muse_model = _circular_moffat_profile(
        x, y, pixel_area,
        muse_results.best_values['dx'], muse_results.best_values['dy'],
        muse_results.best_values['bg'], muse_results.best_values['flux'],
        muse_results.best_values['fwhm'], muse_results.best_values['beta'])

    # Calculate the RMS of the residual image pixels.

    muse_rms_error = np.sqrt(np.mean((muse_values - muse_model)**2))

    # Get the values of the HST pixels at each element of x and y.

    hst_values = hst.data.filled(0.0).ravel()[sorted_indexes]

    # Also attempt to fit a 2D Moffat profile to the selected pixels
    # of the HST image. Again we set the initial guesses for dx and dy
    # to the flux-weighted means of x and y, and the initial guess at
    # the moffat flux to the sum of the pixel values.

    fitmod = Model(_circular_moffat_profile,
                   independent_vars=["x", "y", "pixel_area"])
    fitmod.set_param_hint('fwhm', value=0.5, min=0.0)
    fitmod.set_param_hint('beta', value=10.0, min=0.0)
    fitmod.set_param_hint('dx', value=np.average(x, weights=hst_values))
    fitmod.set_param_hint('dy', value=np.average(y, weights=hst_values))
    fitmod.set_param_hint('flux', value=hst_values.sum(), min=0.0)
    fitmod.set_param_hint('bg', value=0.0)
    fitmod.make_params()

    # Fit the model to the HST data.

    hst_results = fitmod.fit(hst_values, x=x, y=y, pixel_area=pixel_area,
                             weights=10.0)

    # If a hardcopy format has been specified, construct the filename
    # for the saved plot.

    if hardcopy is None or hardcopy == "" or hardcopy == "none":
        plotfile = None
    else:
        prefix = muse.filename.replace(".fits", "")
        plotfile = prefix + "_star_fit." + hardcopy

    # If requested plot the results.

    if display or plotfile is not None:
        _plot_fitted_star_results(muse, hst, radius, muse_results, hst_results,
                                  x, y, muse_values, hst_values,
                                  fig, display, plotfile, nowait, title)

    # Encapsulate the results for return.

    imfit = FittedStarPhotometry(muse, muse_results, hst_results,
                                 muse_rms_error, ra, dec)

    # Write a corrected version of the MUSE image?

    if apply:
        _write_corrected_image(muse, imfit, resample)

    # Return the fitted parameters.

    return imfit


def _write_corrected_image(muse, imfit, resample):
    """Derive corrections from the fitted position errors and
    calibration errors, apply these to the MUSE image, and write
    the resulting image to a FITS file. The name of the output file
    is based on the name of the input file, by replacing its
    ".fits" extension with "_aligned.fits".  If the input muse
    image was not read from a file, then a file called
    "image_aligned.fits" is written in the current directory.

    muse : `mpdaf.obj.Image`
       The image to be corrected.
    imfit : `imphot.FittedPhotometry`
       Fitted image properties.
    resample : False
       If resample=False, then the coordinate reference pixel (CRPIX1,
       CRPIX2) is adjusted to correct the coordinates of the pixels
       without changing any pixel values. If resample=True, then the
       pixel values are resampled to shift the image without changing
       the coordinates of the pixels.
    """

    # Get a corrected version of the MUSE image.

    im = apply_corrections(muse, imfit, resample=resample)

    # Get a prefix for the output filename.

    if muse.filename is None:
        prefix = "muse"
    else:
        prefix = muse.filename.replace(".fits", "")

    # Write the corrected image to disk.

    im.write(prefix + "_aligned.fits")

# Define the class that holds information returned by
# fit_star_photometry().


class FittedStarPhotometry(FittedPhotometry):
    """A class for returning the results of `fit_star_photometry()`.

    Note that both this class and `FittedImagePhotometry` are derived
    from the `FittedPhotometry` class, which contains the essential
    features of the fitted parameters, without features that are
    specific to the fitting method.

    Parameters
    ----------
    muse : `mpdaf.obj.Image`
       The MUSE image that the fit was performed on.
    muse_results : `lmfit.ModelResult`
       The results of fitting a 2D Moffat to a star in the MUSE image.
    hst_results : `lmfit.ModelResult`
       The results of fitting a 2D Moffat to a star in the HST image.
    ra  : float
       The right-ascension used for the origin of the fit (degrees)
    dec : float
       The declination used for the origin of the fit (degrees)
    rms_error : float
       The root-mean square of the pixel residuals of the fit to the
       star in the MUSE image, in the same units as the pixels of the
       original MUSE image.

    Attributes
    ----------
    name : str
       The basename of the MUSE FITS without the .fits extension.
    fit_report : str
       A multi-line string listing verbose reports on both the\n
       MUSE and the HST star fits.
    scale  : `FittedValue`
       The best-fit value and error of the calibration scale\n
       factor, MUSE/HST. This is derived from muse_flux/hst_flux.
    bg  : `FittedValue`
       The best-fit value and error of the calibration offset,\n
       MUSE-HST. This is derived from muse_bg - hst_bg.
    dx : `FittedValue`
       The best-fit value and error of the x-axis pointing offset,\n
       muse_dx-hst_dx (arcsec). This is the distance that features\n
       in the HST image had to be moved to the right, in the\n
       direction of increasing x-axis pixel index, to line them up\n
       with the same features in the MUSE image. One way to\n
       correct the pointing error of the MUSE observation, is to\n
       divide 'dx' by the pixel size along the x-axis (usually\n
       0.2 arcsec), then add the resulting pixel offset to the\n
       CRPIX1 header parameter.
    dy : `FittedValue`
       The best-fit value and error of the y-axis pointing offset,\n
       muse_y-hst_y (arcsec). This is the distance that features\n
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
       PSF (arcsec). This is the same as muse_fwhm.
    beta : `FittedValue`
       The best-fit value and error of the beta parameter of the\n
       Moffat PSF. This is the same as muse_beta.
    rms_error : float
       The root-mean square of the pixel residuals of the fit to\n
       the star in the MUSE image, in the same units as the pixels\n
       of the original MUSE image.

    ra  : float
       The right-ascension used for the origin of the fit (degrees)
    dec : float
       The declination used for the origin of the fit (degrees)

    muse_fwhm : `FittedValue`
       The best-fit value and error of the FWHM of the Moffat\n
       PSF fitted to the star in the MUSE image (arcsec)
    muse_beta : `FittedValue`
       The best-fit value and error of the beta parameter of the\n
       Moffat PSF fitted to the star in the MUSE image.
    muse_dx : `FittedValue`
       The x-axis image offset of the center of the star fitted\n
       in the MUSE image to the coordinate self.ra,self.dec\n
       (arcsec).
    muse_dy : `FittedValue`
       The y-axis image offset of the center of the star fitted\n
       in the MUSE image to the coordinate self.ra,self.dec\n
       (arcsec).
    muse_bg : `FittedValue`
       The fitted flux zero-offset under the fitted Moffat PSF in\n
       the MUSE image of the star.
    muse_flux : `FittedValue`
       The fitted total flux of the Moffat PSF that was fitted to\n
       the star in the MUSE image, using the flux units of the\n
       MUSE image.
    muse.rchi : float
       The reduced chi-squared of the fit of the Moffat PSF to the\n
       star in the MUSE image.

    hst_fwhm : `FittedValue`
       The best-fit value and error of the FWHM of the Moffat\n
       PSF fitted to the star in the HST image (arcsec)
    hst_beta : `FittedValue`
       The best-fit value and error of the beta parameter of the\n
       Moffat PSF fitted to the star in the HST image.
    hst_dx : `FittedValue`
       The x-axis image offset of the center of the star fitted\n
       in the HST image to the coordinate self.ra,self.dec\n
       (arcsec).
    hst_dy : `FittedValue`
       The y-axis image offset of the center of the star fitted\n
       in the HST image to the coordinate self.ra,self.dec\n
       (arcsec).
    hst_bg : `FittedValue`
       The fitted flux zero-offset under the fitted Moffat PSF in\n
       the HST image of the star.
    hst_flux : `FittedValue`
       The fitted total flux of the Moffat PSF fitted to the star\n
       in the HST image, using the flux units of the MUSE image.
    hst_rchi : float
       The reduced chi-squared of the fit of the Moffat PSF to the\n
       star in the HST image.

    """

    def __init__(self, muse, muse_results, hst_results, rms_error, ra, dec):

        # Record the reference RA and Dec of the fits.

        self.ra = ra
        self.dec = dec

        # Record the RMS of the residuals in the MUSE image.

        self.rms_error = rms_error

        # Get the best-fit values for the star in the MUSE image.

        self.muse_fwhm = FittedValue(muse_results.params['fwhm'])
        self.muse_beta = FittedValue(muse_results.params['beta'])
        self.muse_dx = FittedValue(muse_results.params['dx'])
        self.muse_dy = FittedValue(muse_results.params['dy'])
        self.muse_bg = FittedValue(muse_results.params['bg'])
        self.muse_flux = FittedValue(muse_results.params['flux'])
        self.muse_rchi = muse_results.redchi

        # Get the best-fit values for the star in the HST image.

        self.hst_fwhm = FittedValue(hst_results.params['fwhm'])
        self.hst_beta = FittedValue(hst_results.params['beta'])
        self.hst_dx = FittedValue(hst_results.params['dx'])
        self.hst_dy = FittedValue(hst_results.params['dy'])
        self.hst_bg = FittedValue(hst_results.params['bg'])
        self.hst_flux = FittedValue(hst_results.params['flux'])
        self.hst_rchi = hst_results.redchi

        # Calculate the apparent scaling error muse/hst.

        scale = FittedValue(
            value=self.muse_flux.value / self.hst_flux.value,
            stdev=np.sqrt(((self.hst_flux.value * self.muse_flux.stdev)**2 +
                           (self.muse_flux.value * self.hst_flux.stdev)**2)
                          / self.hst_flux.value**4),
            fixed=False)

        # Calculate the apparent background offset error muse-hst.

        bg = FittedValue(
            value=self.muse_bg.value - self.hst_bg.value,
            stdev=np.sqrt(self.muse_bg.stdev**2 + self.hst_bg.stdev**2),
            fixed=False)

        # Calculate the apparent point errors, muse.x - hst.x,
        # and muse.y - hst.y.

        dx = FittedValue(
            value=self.muse_dx.value - self.hst_dx.value,
            stdev=np.sqrt(self.muse_dx.stdev**2 + self.hst_dx.stdev**2),
            fixed=False)
        dy = FittedValue(
            value=self.muse_dy.value - self.hst_dy.value,
            stdev=np.sqrt(self.muse_dy.stdev**2 + self.hst_dy.stdev**2),
            fixed=False)

        # Concatenate the reports from the MUSE and HST fits.

        report = ("The fit of the star in the MUSE image:\n" +
                  muse_results.fit_report() +
                  "\nThe fit of the star in the HST image:\n" +
                  hst_results.fit_report())

        # Set the values of the superclass.

        FittedPhotometry.__init__(self, method="stars", muse=muse,
                                  fit_report=report,
                                  scale=scale, bg=bg, dx=dx, dy=dy,
                                  fwhm=self.muse_fwhm, beta=self.muse_beta,
                                  rms_error=rms_error)


def _plot_fitted_star_results(muse, hst, radius, muse_results, hst_results,
                              x, y, muse_values, hst_values,
                              fig=None, display=True,
                              plotfile=None, nowait=False, title=None):
    """This is a private function of fit_star_photometry() which
    generates a single-page plot of the star fitting results.

    The plot contains 4 graphs. The left-most two plots show the
    MUSE star fit and its residuals, while the right-most pair of
    plots show the HST star fit and its residuals.

    Parameters
    ----------
    muse : `mpdaf.obj.Image`
       The MUSE image.
    hst : `mpdaf.obj.Image`
       The HST image.
    radius : float
       The radius of the area within which the fit was performed.
    muse_results : `lmfit.ModelResult`
       The results of fitting a 2D Moffat to a star in the MUSE image.
    hst_results : `lmfit.ModelResult`
       The results of fitting a 2D Moffat to a star in the HST image.
    x : numpy.ndarray
       The x-axis coordinates the fitted pixels, in increasing order
       of pixel radius from the ra, dec reference position.
    y : numpy.ndarray
       The y-axis coordinates the fitted pixels, in increasing order
       of pixel radius from the ra, dec reference position.
    muse_values : numpy.ndarray
       The values of the MUSE pixels at x, y.
    hst_values : numpy.ndarray
       The values of the HST pixels at x, y.
    fig : matplotlib.figure or None
       The figure in which to plot. If this is None, then
       the default figure, plt.gcf(), is substituted.
    display : bool
       If True (the default), display the plot.
    plotfile : str or None
       An optional filename to plot to, or None. This should have
       a filename suffix that the matplotlib backend understands.
    nowait : bool
       When this argument is False, wait for the user to dismiss
       the plot before returning. This allows the user to interact
       with the plot. When this argument is True, the plot will
       dissapear as soon as another plot is drawn and the user will
       not be able to interact with it.
    title : str or None
       A specific plot title, None or "none" to request the default title,
       or "" if no title is wanted.

    """

    # If no figure has been specified, get the default figure.

    if fig is None:
        fig = plt.gcf()

    # Buffer subsequent graphics until ready to display them.

    plt.ioff()

    # Clear the figure.

    fig.clf()

    # Substitute a default title?

    if title is None or title == "none":
        title = "MUSE image: %s " % basename(muse.filename).replace(".fits", "")

    # Display a plot title?

    if title != "":
        title_y = 0.98 if "\n" in title else 0.95
        fig.suptitle(title, ha="left", x=0.12, y=title_y, fontsize=14)

    # Create a plot grid with 2 rows and 2 columns.

    gs = gridspec.GridSpec(2, 2)

    # Calculate the area of a pixel.

    dy, dx = muse.wcs.get_step(unit=u.arcsec)
    pixel_area = dx * dy

    # Decide the maximum radius to plot.

    xmax = 1.01 * radius

    # Decide where to start the radius axis, giving a small amount
    # of space left of 0 radius, so that there is room to plot the
    # point at zero radius.

    xmin = -xmax / 100.0

    # Create a finely spaced radius axis for plotting the
    # "continuous" model.

    moffat_r = np.linspace(0.0, xmax, 200.0)

    for values, results, col, name in [
            (muse_values, muse_results, 0, "MUSE"),
            (hst_values, hst_results, 1, "HST")]:

        # Calculated the fitted star model at x,y.

        model = _circular_moffat_profile(
            x, y, pixel_area,
            results.best_values['dx'], results.best_values['dy'],
            results.best_values['bg'], results.best_values['flux'],
            results.best_values['fwhm'], results.best_values['beta'])

        # Calculate the model radii of x,y for the MUSE and HST models.

        r = np.sqrt((y - results.best_values['dy'])**2 +
                    (x - results.best_values['dx'])**2)

        # Also calculate the model at the finely sampled radii.

        moffat_values = _circular_moffat_profile(
            moffat_r, 0.0, pixel_area, 0.0, 0.0,
            results.best_values['bg'], results.best_values['flux'],
            results.best_values['fwhm'], results.best_values['beta'])

        # Decide the range of flux-densities to plot.

        ymin = min(0.0, values.min())
        ymax = max(values.max(), moffat_values[0])

        # Render the flux units of the pixels as a latex string,
        # excluding any scale factor.

        units_label = u.Unit(muse.unit / muse.unit.scale).to_string("latex_inline")

        # Matplotlib doesn't yet support \mathring, so replace
        # the anstrom symbol, \mathring{A}, by the smaller
        # unicode angstrom character.

        units_label = units_label.replace(r'\mathring{A}', u"\u00c5")

        # Also remove the excessive '\,' spaces that astropy.units
        # adds after exponents.

        units_label = re.sub(r'(\^\{-?[0-9]+\})\\,', r'\1', units_label)

        # To ensure that the Y-axis flux and residual flux labels
        # don't end up dominated by lots of zeros before the decimal
        # point, calculate the power-of-10 scale factor needed to
        # scale the maximum flux value to have 2 digits before the
        # decimal point.

        y_multiplier = 10.0**(-np.floor(np.log10(abs(ymax))) + 1)

        # Compose a units label for labelling the Y-axis.

        y_units_label = "%s%s" % (u.Unit(muse.unit.scale / y_multiplier).to_string("latex_inline"), units_label)

        # Rescale the y-axis values and the y-axis range.

        ymin *= y_multiplier
        ymax *= y_multiplier
        values *= y_multiplier
        moffat_values *= y_multiplier
        model *= y_multiplier

        # Plot the actual PSF and the model Moffat function.

        ax = fig.add_subplot(gs[0, col])
        ax.set_autoscale_on(False)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        if col == 0:
            ax.set_ylabel(y_units_label)
        traces = []; labels = []
        ax.plot(r, values, ls="", marker=".")
        ax.plot(moffat_r, moffat_values)

        # Label the contents of the graph.

        ax.text(0.2, 0.95, (("%s pixels plotted over a\n" % name) +
                            "Moffat model of total flux:"),
                transform=ax.transAxes, ha="left", va="top")
        ax.text(0.28, 0.82,
                u"%.3g %s\nfwhm: %.2g arcsec\nMoffat beta: %.2g" % (
                    results.best_values['flux'] * muse.unit.scale, units_label,
                    results.best_values['fwhm'], results.best_values['beta']),
                transform=ax.transAxes, ha="left", va="top")

        # Plot the residuals between the actual PSF and the Moffat PSF.

        ax = fig.add_subplot(gs[1, col])
        ax.set_autoscale_on(False)
        ax.set_xlim(xmin, xmax)
        if col == 0:
            ax.set_ylabel(y_units_label)
        ax.set_xlabel("Radius (arcsec)")
        residuals = (values - model)
        res_min = residuals.min()
        res_max = residuals.max()
        res_margin = (res_max - res_min) * 0.5
        ax.set_ylim(res_min - res_margin, res_max + res_margin)
        ax.plot(r, residuals, ls="", marker=".")
        ax.text(0.5, 0.95, "%s residuals (pixels - model)" % name,
                transform=ax.transAxes, ha="center", va="top")

    # Display the plot?

    if display:
        if not nowait:
            plt.show()
        else:
            plt.ion()
            plt.pause(0.5)

    # Record the plot to a file?

    if plotfile is not None:
        fig.savefig(plotfile, orientation="landscape", dpi=600)


def _circular_moffat_profile(x, y, pixel_area, dx, dy, bg, flux, fwhm, beta):
    """Return the flux/pixel of a circularly symmetric Moffat profile at a
    specified x, y position.

    Parameters
    ----------
    x    :  float
       The X-axis position at which to evaluate the moffat
       function (arcsec).
    y    :  float
       The Y-axis position at which to evaluate the moffat
       function (arcsec).
    pixel_area : float
       The area of the pixel being sampled (arcsec**2).
       The returned flux density will be the brightness at
       x,y scaled by this area.
    dx     :  float
       The X-axis offset of the center of the Moffat function
       from the origin of X.
    dy     :  float
       The Y-axis offset of the center of the Moffat function
       from the origin of Y.
    bg     :  float
       The flux offset to which the Moffat function is to be added.
    flux   :  float
       The integrated flux under the Moffat function.
    fwhm   :  float
       The full-width at half-maximum of the Moffat function (arcsec).
    beta   : float
       The term due to scattering in the atmosphere that widens
       the wings of the PSF compared to a Gaussian. Values above
       about 5.0 make the Moffat function very similar to a
       Gaussian profile. Values below this make the profile fall-off
       more gradually than a gaussian, particularly outside the FWHM
       of the profile.

    Returns
    -------
    out : numpy.ndarray
       The value of the Moffat profile at x,y, in units of flux/pixel.

    """

    # A 2D Moffat function is defined as follows:
    #
    #   y(x,y) = flux * (beta-1)/(pi*a**2) / (1 + (x**2 + y**2) / a**2)**beta
    #
    # Calculate a**2 from the FWHM and beta.

    asq = fwhm**2 / 4.0 / (2.0**(1.0 / beta) - 1.0)

    # Compute the peak value.

    peak = flux * (beta - 1.0) / (np.pi * asq)

    # Compute the brightness of the Moffat function at x,y,
    # scaled by the specified pixel area, to convert the brightness
    # to a flux density per pixel.

    return bg + pixel_area * peak / (1.0 + ((x - dx)**2 + (y - dy)**2) / asq)**beta


class FitStarPhotometryMP(_FitPhotometryMP):
    """A multiprocessing iterator that creates a pool or worker
    processes to repeatedly call `fit_star_photometry()` for each
    of a list of MUSE files, returning the results via the
    iterator as they become available.

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
    nworker : int
       The number of worker processes to use. The default is
       0, which creates multiprocessing.cpu_count() processes.
       Alternatively, if a negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.
    kwargs : dict
       A dictionary of keyword/value arguments to be passed to
       fit_star_photometry(). This should always include a value
       for the star=(ra,dec,radius) argument.

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

        super(FitStarPhotometryMP, self).__init__(
            hst_filename, muse_filenames,
            cmd_fn=fit_star_photometry, cmd_kwargs=kwargs, nworker=nworker)

    def next(self):
        """Return the results from the next image in the list of input MUSE
        images.

        Returns
        -------
        out : `FittedStarPhotometry`
           The fitting results from the next image.

        """

        return super(FitStarPhotometryMP, self).next()
