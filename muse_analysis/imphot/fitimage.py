from __future__ import division
from os.path import basename
import numpy as np
from numpy import ma
import astropy.units as u
from astropy.io import fits
from scipy.ndimage.morphology import (grey_opening, grey_dilation, binary_erosion)
from scipy.ndimage.filters import convolve
from lmfit import Model
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from mpdaf.obj import Image

from . import ds9regions
from .core import (UserError, FittedValue, FittedPhotometry, HstFilterInfo,
                   image_grids_aligned, _default_hst_fwhm, _default_hst_beta)
from .mp import _FitPhotometryMP

__all__ = ['fit_image_photometry', 'FittedImagePhotometry',
           'FitImagePhotometryMP']

# Set the width of the chamfer to apply to image masks. This should be
# an odd valued integer, to make it symmetric about a central pixel.

_bevel_width = 9   # Integer pixels

def fit_image_photometry(hst, muse, regions=None, fix_scale=None,
                         fix_bg=None, fix_dx=None, fix_dy=None,
                         fix_fwhm=None, fix_beta=None,
                         hst_fwhm=_default_hst_fwhm,
                         hst_beta=_default_hst_beta,
                         margin=2.0, segment=False, display=False,
                         nowait=False, hardcopy=None, title=None,
                         star=None, save=False, fig=None):

    """Given a MUSE image and an HST image that has been regridded and aligned
    onto the same coordinate grid as the MUSE image, use the HST image as a
    calibrated reference to fit for the flux-scale, the FWHM and beta
    parameter of a Moffat PSF, the position offset, and the background
    offset of the MUSE image relative to the HST image.

    More specifically, it uses the HST image to reproduce the MUSE
    image by convolving the HST image with a PSF that best reproduces
    the resolution of the MUSE image, scales and offsets the HST image
    with numbers that best reproduce the fluxes of the MUSE image, and
    shifts the HST image by an amount that best lines up features in
    the two images. The PSF parameters, calibration factors and
    pointing offsets that best reproduce the MUSE image, are the
    outputs of this algorithm.

    Parameters
    ----------
    hst : `mpdaf.obj.Image` or filename
       An HST image of the same area as the MUSE image, and that has been
       gridded onto the same pixel grid as the MUSE image. This can be
       given as an MPDAF Image object, or by the filename of the FITS file.
    muse : `mpdaf.obj.Image` or filename
       The MUSE image to be characterized. This can be
       given as an MPDAF Image object, or by the filename of the FITS file.
    regions : None, str, or iterable
       This can be None, to indicate that no regions are needed,
       "star" to restrict the fit to pixels within any region defined
       by the star=(ra,dec,radius) argument, "notstar" to restrict the
       fit to pixels outside any region defined the by
       star=(ra,dec,radius) argument, the name of a filename of the
       ds9 region file, or an iterable that returns successive lines
       of a ds9 region file.

       The regions can be used to exclude problematic areas of an
       image or sources that would degrade the global PSF fit, such as
       saturated stars, stars with significant proper motion, and
       variable sources.  Alternatively it can be used to restrict the
       fit to one or more objects by masking everything except small
       regions around these objects.

       Only ds9 circle, ellipse and box regions are supported. Other
       region types and most of the configuration parameters found in
       ds9 region files, are simply ignored. If the file is created
       within ds9, be careful to tell ds9 whether you want the region
       to be included or excluded. Also be careful to request either
       fk5 or physical pixel coordinates, because other coordinate
       systems are not supported.
    fix_scale : float or None
       The calibration scale factor, (MUSE_flux / HST_flux) is fixed
       to the specified value while fitting, unless the value is None.
    fix_bg : float or None
       The calibration zero-offset, (MUSE_flux - HST_flux) is fixed
       to the specified value while fitting, unless the value is None.
    fix_dx : float or None
       The x-axis pointing offset, (MUSE_x - HST_x) is fixed
       to the specified value while fitting, unless the value is None.
    fix_dy : float or None
       The y-axis pointing offset, (MUSE_y - HST_y) is fixed
       to the specified value while fitting, unless the value is None.
    fix_fwhm : float or None
       The FWHM of the Moffat PSF is fixed to the specified value
       while fitting, unless the value is None.
    fix_beta : float or None
       The beta exponent of the Moffat PSF is fixed to the specified value
       while fitting, unless the value is None.
    hst_fwhm : float
       The FWHM of a Moffat model of the effective PSF of the HST.
       The default value that is used if this parameter is not
       specified, came from Moffat fits to stars within HST UDF
       images. To obtain the closest estimate to the dithered
       instrumental PSF, these fits were made to images with the
       smallest available pixel size (30mas).
    hst_beta : float
       The beta parameter of a Moffat model of the effective PSF of
       the HST.  The default value came from Moffat fits to stars
       within HST UDF images, as described above for the hst_fwhm
       parameter.
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
       "jpg", "png" or "eps". Plots of the image fit will be written to
       a filename that starts with the name of the MUSE input file,
       after removing any .fits suffix, followed by "_image_fit.<suffix>".
    nowait : bool
       When this argument is False, wait for the user to dismiss
       the plot before returning. This allows the user to interact
       with the plot. When this argument is True, the plot will
       dissapear as soon as another plot is drawn and the user will
       not be able to interact with it.
    title : str or None
       A specific plot title, or None to request the default title.
       Specify "" if no title is wanted.
    star : (float, float, float)
       This option can be used to restrict the fitting procedure to
       a single star within the MUSE image. Its arguments are the
       right ascension and declination of the star in decimal degrees,
       and the radius of the area over which to perform the fit, in
       arcseconds. Beware that this option overrides the
       regions argument.
    save : bool
       If True, save the result images of each input image to
       FITS files.

    Returns
    -------
    out : `FittedImagePhotometry`
       An object that contains the fitted parameters and a textual
       report about the fit.

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

    # Get the Y-axis and X-axis dimensions of the MUSE pixels in arcseconds.

    dy, dx = muse.get_step(unit=u.arcsec)

    # Get the masks of the MUSE and HST images.

    if muse.data.mask is not ma.nomask:
        muse_mask = muse.data.mask
    else:
        muse_mask = ~np.isfinite(muse.data.data)

    if hst.data.mask is not ma.nomask:
        hst_mask = hst.data.mask
    else:
        hst_mask = ~np.isfinite(hst.data.data)

    # Before comparing the images we need to make sure that any areas
    # that are masked in one image are also masked in the other, so
    # obtain the union of the MUSE and HST masks.

    mask = np.logical_or(muse_mask, hst_mask)
    del muse_mask
    del hst_mask

    # If requested, mask areas that don't contain significant objects
    # by masking pixels that are below the median value in a
    # morphologically opened version of the HST image.

    if segment:
        tmp = hst.data.filled(0.0)
        kernel = np.array([[0,1,0],
                           [1,1,1],
                           [0,1,0]])
        tmp = grey_opening(tmp, structure=kernel)
        tmp = grey_dilation(tmp, structure=kernel)
        tmpmask = tmp < np.median(tmp)
        tmpmask = binary_erosion(tmpmask, structure=kernel)
        mask = np.logical_or(mask, tmpmask)

    # If we have been passed a circular region around a star,
    # extract its center coordinate and radius.

    if star is not None:
        try:
            ra, dec, radius = star
        except Exception:
            raise UserError("The star argument should be a tuple of (ra,dec,radius)")
        if radius <= 0.0:
            raise UserError("The radius given to the star argument must be > 0")

    # Have we been asked to restrict the fit to pixels inside the
    # circular star region (if any)?

    if regions == "star":
        if star is not None:
            regions = ["fk5; circle(%g, %g, %g\")" % (ra, dec, radius)]
        else:
            regions = None

    # Have we been asked to restrict the fit to pixels outside the
    # circular star region (if any)?

    elif regions == "notstar":
        if star is not None:
            regions = ["fk5; -circle(%g, %g, %g\")" % (ra, dec, radius)]
        else:
            regions = None

    # If the user wants any regions masked, mask them. The masking
    # code needs the mask to be part of an MPDAF Image, so temporarily
    # replace the HST mask with the above mask.

    if regions is not None:
        try:
            original_mask = hst.mask
            hst.mask = mask
            _mask_ds9regions(hst, regions)
            mask = hst.mask
        finally:    # Always reinstall the original HST mask.
            hst.mask = original_mask

    # Copy both image arrays into masked array containers that use the
    # above mask.

    mdata = ma.array(data=muse.data.data, mask=mask, copy=False)
    hdata = ma.array(data=hst.data.data, mask=mask, copy=False)

    # Create boolean arrays that indicate which rows and columns
    # contain no unmasked values.

    masked_rows = np.apply_over_axes(np.logical_and.reduce, mask, 1).ravel()
    masked_cols = np.apply_over_axes(np.logical_and.reduce, mask, 0).ravel()

    # Get the indexes of all rows and columns that contain at least one
    # unmasked element.

    used_rows = np.where(~masked_rows)[0]
    used_cols = np.where(~masked_cols)[0]

    # Create a list of slices to select the above area that contains
    # unmasked pixels.

    crop_indexes = [slice(used_rows.min(), used_rows.max()-1),
                    slice(used_cols.min(), used_cols.max()-1)]

    # Crop the HST and MUSE images and the mask.

    mdata = mdata[crop_indexes]
    hdata = hdata[crop_indexes]
    mask = mask[crop_indexes]

    # Get the median of the unmasked parts of the MUSE image, to use
    # as an initial estimate of the constant part of the background.

    subtracted = np.ma.median(mdata)

    # Get ndarray versions of the above arrays with masked pixels
    # filled with zeros. Note that the choice of zeros (as opposed
    # to median values) prevents these pixels from biasing the
    # fitted values of the image offset and scale factors.
    #
    # Before doing this, subtract the estimated MUSE background from
    # the MUSE image, such that when the masked areas are replaced
    # with zeros, they aren't significantly different from the normal
    # background of the image.

    mdata = ma.filled(mdata - subtracted, 0.0)
    hdata = ma.filled(hdata, 0.0)

    # When a position offset in the image plane is performed by
    # applying a linearly increasing phase offset in the Fourier
    # domain of an FFT, the result is that features in the image
    # that get shifted off the right edge of the image reappear at the
    # left edge of the image, and vice versa. By appending a margin of
    # zeros to the right edge, and later discarding this margin, we
    # ensure that features that are shifted off the left or right
    # edges of the image end up in this discarded area, rather than
    # appearing at the opposite edge of the visible image. Appending a
    # similar margin to the Y axis has the same effect for vertical
    # shifts.

    shape = np.asarray(mdata.shape) + np.ceil(np.abs(
        margin / muse.get_step(unit=u.arcsec) ) ).astype(int)

    # Round the image dimensions up to integer powers of two, to
    # ensure that an efficient FFT implementation is used.

    shape = (2**np.ceil(np.log(shape)/np.log(2.0))).astype(int)

    # Extract the dimensions of the expanded Y and X axes.

    ny,nx = shape

    # Compute the slice needed to extract the original area from
    # expanded arrays of the above shape.

    sky_slice = [slice(0,mdata.shape[0]), slice(0,mdata.shape[1])]

    # Zero-pad the MUSE image array to have the new shape.

    tmp = np.zeros(shape)
    tmp[sky_slice] = mdata
    mdata = tmp

    # Zero-pad the HST image array to have the new shape.

    tmp = np.zeros(shape)
    tmp[sky_slice] = hdata
    hdata = tmp

    # Pad the mask array to have the same dimensions, with padded
    # elements being masked.

    tmp = np.ones(shape, dtype=bool)
    tmp[sky_slice] = mask
    mask = tmp

    # Compute a multiplicative pixel scaling image that inverts the
    # mask to be True for pixels that are to be kept, and zero for
    # pixels that are to be removed, then smoothly bevel the edges of
    # the areas of ones to reduce the effects of sharp edges on the
    # fitted position offsets.

    weight_img = _bevel_mask(~mask, _bevel_width)

    # Also obtain the FFT of the mask for fitting the background flux
    # offset in the Fourier plane.

    weight_fft = np.fft.rfft2(weight_img)

    # Scale the MUSE and HST images by the weighting image.

    mdata *= weight_img
    hdata *= weight_img

    # Obtain the Fourier transforms of the MUSE and HST images. The
    # Fourier transform is hermitian, because the images have no
    # imaginary parts, so we only calculate half of the Fourier
    # transform plane by using rfft2 instead of fft2().

    hfft = np.fft.rfft2(hdata)
    mfft = np.fft.rfft2(mdata)

    # Scale the FFT of the MUSE image by the blackman
    # anti-aliasing window function which is assumed to have been
    # applied by Image.resample() to decimate the HST image to the
    # MUSE image resolution.

    _apply_resampling_window(mfft, shape)

    # Convolve the muse function with the PSF of the original HST
    # image.

    _convolve_hst_psf(mfft, shape, dy, dx, hst_fwhm, hst_beta)

    # Compute the spatial frequencies along the X and Y axes
    # of the above power spectra.

    fx = np.fft.rfftfreq(nx, dx)
    fy = np.fft.fftfreq(ny, dy)

    # Get 2D grids of the x,y spatial-frequency coordinates of each pixel
    # in the power spectra.

    fx2d,fy2d = np.meshgrid(fx,fy)

    # The initial guess for the factor to multiply the HST fluxes
    # by, is the ratio of the values at the origins of the MUSE
    # and HST power spectra. The values at the origin represent
    # the sums of the flux in the MUSE and HST images.

    scale_guess = np.abs(mfft[0,0]) / np.abs(hfft[0,0])

    # Calculate the frequency interval of the FFTs along the
    # X and Y axes.

    dfx = 1.0 / (nx * dx)
    dfy = 1.0 / (ny * dy)

    # Get a 2D array of the radius of each image pixel center relative
    # to pixel [0,0] (the spatial origin of the FFT algorithm).

    rsq = np.fft.fftfreq(nx, dfx)**2 + \
          np.fft.fftfreq(ny, dfy)[np.newaxis,:].T**2

    # Wrap the modeling function in an object for the least-squares fitting
    # algorithm.

    fitmod = Model(_xy_moffat_model_fn,
                   independent_vars=["fx", "fy", "rsq", "hstfft", "wfft",
                                     "subtracted"])

    # Describe each of the parameters of the model to the least-squares fitter.

    if fix_bg is None:
        fitmod.set_param_hint('bg', value=subtracted)
    else:
        fitmod.set_param_hint('bg', value=fix_bg,
                                vary=False)
    if fix_scale is None:
        fitmod.set_param_hint('scale', value=scale_guess, min=0.0)
    else:
        fitmod.set_param_hint('scale', value=fix_scale,
                                min=0.0, vary=False)
    if fix_fwhm is None:
        fitmod.set_param_hint('fwhm', value=0.5, min=0.0)
    else:
        fitmod.set_param_hint('fwhm', value=fix_fwhm,
                                min=0.0, vary=False)
    if fix_dx is None:
        fitmod.set_param_hint('dx', value=0.0)
    else:
        fitmod.set_param_hint('dx', value=fix_dx, vary=False)

    if fix_dy is None:
        fitmod.set_param_hint('dy', value=0.0)
    else:
        fitmod.set_param_hint('dy', value=fix_dy, vary=False)

    if fix_beta is None:
        fitmod.set_param_hint('beta', value=5.0)
    else:
        fitmod.set_param_hint('beta', value=fix_beta, vary=False)
    fitmod.make_params()

    # Fit for any parameters that have been allowed to vary.

    results = fitmod.fit(mfft.ravel().view(dtype=float),
                         fx=fx2d.ravel(), fy=fy2d.ravel(),
                         rsq=rsq, hstfft=hfft.ravel(),
                         wfft=weight_fft.ravel(), subtracted=subtracted)

    # Compute the FFT of the modified HST image.

    hfft = _xy_moffat_model_fn(fx2d, fy2d, rsq, hfft, weight_fft, subtracted,
                               results.best_values['dx'],
                               results.best_values['dy'],
                               results.best_values['bg'],
                               results.best_values['scale'],
                               results.best_values['fwhm'],
                               results.best_values['beta']).view(dtype=complex)

    # Invert the best-fit FFTs to obtain the best-fit MUSE and HST images.

    muse_im = np.fft.irfft2(mfft)[sky_slice]
    hst_im = np.fft.irfft2(hfft)[sky_slice]

    # Compute the root-mean square of the residuals.

    rms_error = ma.sqrt((ma.array(muse_im-hst_im, mask=mask[sky_slice])**2).mean())

    # Extract relevant results for return.

    imfit = FittedImagePhotometry(muse, results=results, rms_error=rms_error)

    # If a hardcopy format has been specified, construct the filename
    # for the saved plot.

    if hardcopy is None or hardcopy == "":
        plotfile = None
    else:
        prefix = basename(muse.filename).replace(".fits","")
        plotfile = prefix + "_image_fit." + hardcopy

    # If needed, generate the images needed for plotting and saving.

    if display or plotfile is not None or save:
        (muse_ft, hst_ft) = _generate_fft_images(mfft, hfft)

        # Display the images?

        if display or plotfile is not None:
            _plot_fitted_image_results(muse, imfit, muse_im, muse_ft,
                                       hst_im, hst_ft, fig, display,
                                       plotfile, nowait, title)

        # Save the images?

        if save:
            _save_fitted_images(muse, muse_im, muse_ft,
                                hst, hst_im, hst_ft, crop_indexes)

    # Return the results in an object.

    return imfit

# Define the class that holds information returned by
# fit_image_photometry().

class FittedImagePhotometry(FittedPhotometry):
    """The class of the object that `fit_image_photometry()` returns.

    Note that both this class and `FittedStarPhotometry` are derived
    from the `FittedPhotometry` class, which contains the essential
    features of the fitted parameters, without features that are
    specific to the fitting method.

    Parameters
    ----------
    muse : `mpdaf.obj.Image`
       The MUSE image that the fit was performed on.
    results : `lmfit.ModelResult`
       The model fitting results.
    rms_error : float
       The root-mean square of the residual image pixels, in the same
       units as the pixels of the original MUSE image.

    Attributes
    ----------
    name : str
       The basename of the MUSE FITS without the .fits extension.
    fit_report : str
       A printable report on the fit from the least-squares\n
       fitting function.
    scale  : `FittedValue`
       The best-fit value and error of the calibration scale\n
       factor, (MUSE.flux / HST.flux).
    bg  : `FittedValue`
       The best-fit value and error of the calibration offset,\n
       (MUSE.flux - HST.flux).
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
       The best-fit value and error of the FWHM of the Moffat PSF.
    beta : `FittedValue`
       The best-fit value and error of the beta parameter of the\n
       Moffat PSF
    rchi : float
       The reduced chi-squared value of the fit. Beware that the\n
       absolute value of this number is not very useful, because\n
       the fit is performed to the complex pixels of zero-padded\n
       FFTs, which are not independent observables, rather than to\n
       image pixels.
    rms_error : float
       The root-mean square of the residual image pixels, in the\n
       same units as the pixels of the original MUSE image.

    """

    def __init__(self, muse, results, rms_error):
        FittedPhotometry.__init__(self, method="image", muse=muse,
                                  fit_report=results.fit_report(),
                                  scale=FittedValue(results.params['scale']),
                                  bg=FittedValue(results.params['bg']),
                                  dx=FittedValue(results.params['dx']),
                                  dy=FittedValue(results.params['dy']),
                                  fwhm=FittedValue(results.params['fwhm']),
                                  beta=FittedValue(results.params['beta']),
                                  rms_error=rms_error)
        self.rchi = results.redchi

def _xy_moffat_model_fn(fx, fy, rsq, hstfft, wfft, subtracted, dx, dy,
                       bg, scale, fwhm, beta):

    """This function is designed to be passed to lmfit to fit the FFT of
    an HST image to the FFT of a MUSE image on the same coordinate
    grid.

    It takes the FFT of an HST image and changes it as follows:

    1. It multiplies the FFT by the scale parameter, to change its flux
       calibration.
    2. It multiplies the FFT by the FFT of a circularly symmetric moffat
       function of the specified fwhm and beta parameter. This is
       equivalent to convolving the original HST image with a
       Moffat PSF with the specified characteristics.
    3. It multiplies the FFT by phase gradients that are equivalent to
       shifting the original HST image by dx and dy arcseconds, to
       correct any position offset.
    4. It adds bg*wfft to change the zero-offset of the original HST
       image by bg.

    Parameters:
    -----------
    fx : numpy.ndarray
       The X-axis spatial-frequency coordinate (cycles/arcsec) of each
       pixel in the HST FFT.
    fy : numpy.ndarray
       The Y-axis spatial-frequency coordinate (cycles/arcsec) of each
       pixel in the HST FFT.
    rsq : numpy.ndarray
       The radius-squared of each pixel of an image of the same
       dimensions as the image that was forier transformed to obtain
       hstfft, relative to the center of pixel [0,0], in the same
       spatial units as fwhm. Note that the x offsets of each pixel
       along the x axis, relative to pixel index 0, can be calculated
       using np.fft.fftfreq(nx,1.0/(nx*dx)), where nx is the number of
       pixels along the x axis, and dx is the pixel width in
       arcseconds.  Similarly for the y offsets of each pixel.
    hstfft : numpy.ndarray
       The FFT of the HST image.
    wfft : numpy.ndarray
       The FFT of the weighting array that has been applied to the
       pixels of the image.
    subtracted : float
       A constant background value that has was already subtracted from
       the MUSE image before it was FFT'd.
    dx: float
       The position offset along the X axis of the image.
    dy: float
       The position offset along the Y axis of the image.
    bg : float
       The background offset to add to the HST image (in the Fourier plane
       this is added to the pixel at the origin, which contains the sum
       of the flux in the image plane.)
    scale : float
       The scale factor to multiply the HST FFT pixels by.
    fwhm : float
       The full-width at half maximum of the Moffat function in the image
       plane that will be convolved with the HST image.
    beta : float
       The term due to scattering in the atmosphere that widens
       the wings of the PSF compared to a Gaussian. A common
       choice for this valus is 2.0.

    Returns:
    --------
    out : numpy.ndarray
       The modified value of the HST FFT to be compared to the
       MUSE FFT. The fitter can't cope with complex elements,
       so this function returns an array of two float elements (real, imag)
       for each element of the FFT. To convert this back to a complex array,
       note that you can follow it with .view(dtype=complex)

    """

    # A center-normalized Moffat function is defined as follows:
    #
    #   y(x,y) = 1 / (1 + (x**2 + y**2) / a**2)**beta
    #
    # Calculate a**2 from the FWHM and beta.

    asq = fwhm**2 / 4.0 / (2.0**(1.0 / beta) - 1.0)

    # Compute an image of a Moffat function centered at pixel 0,0.

    im = 1.0 / (1.0 + rsq / asq)**beta

    # Obtain the discrete Fourier Transform of the Moffat function.
    # The function is even, meaning that its Fourier transform is
    # entirely real, so also discard the imaginary parts.

    moffat_ft = np.real(np.fft.rfft2(im))

    # Normalize it to have unit volume in the image plane.

    moffat_ft /= moffat_ft[0,0]

    # If the HST FFT array has been flattened to make it compatible with
    # lmfit, do the same to moffat_ft.

    if hstfft.ndim == 1:
        moffat_ft = moffat_ft.ravel()

    # Precompute the image shifting coefficients.

    argx = -2.0j * np.pi * dx
    argy = -2.0j * np.pi * dy

    # Create the model to compare with the MUSE FFT. This is the HST FFT
    # scaled by the fitted scaling factor, smoothed to a lower resolution
    # by the above 2D Moffat function, and shifted by dx and dy.

    model = (bg - subtracted) * wfft + hstfft * scale * moffat_ft * np.exp(argx*fx + argy*fy)

    # The model-fitting function can't handle complex numbers, so
    # return the complex FFT model as an array of alternating real and
    # imaginary floats.

    return model.view(dtype=float)

def _plot_2d_array(data, axes, vmin=None, vmax=None, pixw=None, pixh=None,
                   cmap=None, xlabel=None, ylabel=None, title=None,
                   title_fs=10, title_color="black", xtick=None, ytick=None,
                   axis_bg="#cccccc"):

    # Substitute defaults for omitted parameters.

    if cmap is None:
        cmap = plt.cm.gray

    if pixw is None:
        pixw = data.shape[1]
    if pixh is None:
        pixh = data.shape[0]

    if vmin is None or vmax is None:
        std = data.std()
        mean = data.mean()
        vmin = mean - 0.4 * std
        vmax = mean + 2.0 * std

    # Get the axis extent of the image.

    im_lft = 0
    im_rgt = data.shape[1]
    im_bot = 0
    im_top = data.shape[0]

    # Axis limits.

    ax_lft = 0
    ax_rgt = pixw
    ax_bot = 0
    ax_top = pixh

    # Create the graph.

    axes.set_axis_bgcolor(axis_bg)
    axes.set_autoscale_on(False)
    axes.set_xlim(ax_lft, ax_rgt)
    axes.set_ylim(ax_bot, ax_top)
    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    if title is not None:
        axes.set_title(title, fontsize=title_fs, color=title_color)
    axes.set_aspect('equal')
    if xtick is not None:
        axes.set_xticks(np.arange(ax_lft, ax_rgt, xtick))
    if ytick is not None:
        axes.set_yticks(np.arange(ax_bot, ax_top, ytick))
    axes.imshow(data, origin="lower",
                vmin=vmin, vmax=vmax,
                extent=(im_lft,im_rgt,im_bot,im_top), cmap=cmap,
                interpolation="nearest")

def _bevel_mask(mask, width):
    """Return a floating point image that is a smoothed version of a
    boolean mask array. It is important that pixels that are zero in
    the input masked array remain zero in the smoothed version. If we
    just convolved the mask with a smoothing kernel, this would smear
    non-zero values into the region of the zeros, so instead we first
    erode the edges of the masked regions using a square kernel, then
    apply a square smoothing kernel of the same size to smooth down to
    the edges of the original masked areas.

    Parameters
    ----------
    mask   :   numpy.ndarray
       A 2D array of bool elements.
    width  :   int
       The width of the bevel in pixels.

    Returns
    -------
    out    :   numpy.ndarray
       A 2D array of float elements.

    """

    # First shave off 'width' pixels from all edges of the mask,
    # and return the result as a floating point array.

    im = binary_erosion(mask, structure=np.ones((width, width))).astype(float)

    # Compute a [width,width] smoothing convolution mask.

    w = np.blackman(width+2)[1:width+1]
    m = w * w[np.newaxis,:].T
    m /= m.sum()

    # Smooth the eroded edges of the mask.

    return convolve(im, m)

def _apply_resampling_window(imfft, shape):

    """Multiply the FFT of an image with the blackman anti-aliasing window
    that would be used by default by the MPDAF Image.resample()
    function for the specified image grid. The multiplication is performed
    in-place, so this function has no return value.

    Parameters
    ----------
    imfft  : numpy.ndarray
       The FFT of the resampled image. This must just be the positive
       frequency half of the FFT of a real image, as returned by
       numpy.fft.rfft2(im).
    shape  : (ny,nx)
       The dimensions of the image that was transformed to yield
       imfft. Note that while ny can be obtained from imfft.shape[0],
       nx can't because it could either be 2*imfft.shape[0]-2 or
       2*imfft.shape[0]-1, depending on whether nx was even or odd.

    """

    # Create an array which, for each pixel in the FFT image, holds
    # the radial spatial-frequency of the pixel center, divided by the
    # Nyquist folding frequency (ie. 0.5 cycles/image_pixel).

    f = np.sqrt((np.fft.rfftfreq(shape[1]) / 0.5)**2 +
                (np.fft.fftfreq(shape[0]) / 0.5)[np.newaxis,:].T**2)

    # Scale the FFT by the window function, calculating the blackman
    # window function as a function of frequency divided by its cutoff
    # frequency.

    imfft *= np.where(f <= 1.0, 0.42 + 0.5 * np.cos(np.pi * f) + 0.08 * np.cos(2*np.pi * f), 0.0)

def _convolve_hst_psf(imfft, shape, dy, dx, fwhm, beta):
    """Convolve an image of MUSE pixel resolution with an HST PSF,
    modelled as a Moffat function. The convolution is performed
    in-place on the imfft argument, so this function does not have a
    formal return value.

    Parameters
    ----------
    imfft  : numpy.ndarray
       The FFT of the resampled image. This must just be the positive
       frequency half of the FFT of a real image, as returned by
       numpy.fft.rfft2(im).
    shape  : (ny,nx)
       The dimensions of the image that was transformed to yield
       imfft. Note that while ny can be obtained from imfft.shape[0],
       nx can't because it could either be 2*imfft.shape[0]-2 or
       2*imfft.shape[0]-1, depending on whether nx was even or odd.
    dy     : The Y-axis resolution of the resampled image (arcsec).
    dx     : The X-axis resolution of the resampled image (arcsec).
    fwhm   : The FWHM of the Moffat model of the HST PSF.
    beta   : The beta parameter of the Moffat model of the HST PSF.

    """

    # A normalized Moffat function as a function of the radius, r,
    # is defined as follows:
    #
    #   moffat(r) = 1 / (1 + r**2 / a**2)**beta
    #
    # Calculate a**2 from the FWHM and beta.

    asq = fwhm**2 / 4.0 / (2.0**(1.0 / beta) - 1.0)

    # Get the dimensions of the MUSE image array.

    ny, nx = shape

    # The FWHM of the HST PSF is badly undersampled at the MUSE
    # resolution, so sample it at an integer fraction of the MUSE
    # pixel size.

    mag = 8

    # Get a 2D array of the radius of each image pixel center
    # relative to the center of pixel image 0,0 (the origin of the
    # FFT algorithm).

    dfx = 1.0 / (nx * dx)
    dfy = 1.0 / (ny * dy)
    mrsq = np.fft.fftfreq(nx * mag, dfx)**2 + \
           np.fft.fftfreq(ny * mag, dfy)[np.newaxis,:].T**2

    # Compute an image of a Moffat function centered at pixel 0,0.

    im = 1.0 / (1.0 + mrsq / asq)**beta

    # Obtain the discrete Fourier Transform of the Moffat function.
    # The function is even, meaning that its Fourier transform is
    # entirely real, so also discard the imaginary parts.

    moffat_ft = np.real(np.fft.rfft2(im))
    del(im)

    # Normalize it to have unit volume in the image plane.

    moffat_ft /= moffat_ft[0,0]

    # Multiply the MUSE FFT by the FFT of the Moffat function.

    imfft[0:ny//2, 0:nx//2] *= moffat_ft[0:ny//2, 0:nx//2]
    imfft[ny-ny//2:ny, 0:nx//2] *= moffat_ft[mag*ny-ny//2:mag*ny, 0:nx//2]
    del(moffat_ft)

def _mask_inside_ds9region(im, reg):
    """A private function of _mask_ds9regions() that masks all pixels
    within a region, regardless of whether that region is marked to
    be included or excluded.

    Parameters
    ----------
    im   : `mpdaf.obj.Image`
       The image to be masked.
    reg  : `ds9regions.Ds9Region`
       The region to be masked.
    """

    # Accomodate the oordinate system that the region described in.

    if reg.system == "fk5":
        units = u.deg
        if hasattr(reg, "pa"):
            pa = reg.pa + im.get_rot()
        center = np.array([reg.y, reg.x])
    elif reg.system == "physical":
        units = None
        if hasattr(reg, "pa"):
            pa = reg.pa
        center = np.array([reg.y, reg.x]) - 1.0
    else:
        return

    # Use the MPDAF function that masks the shape of the region.

    if reg.shape == "circle":
        im.mask_region(center=center, radius=reg.radius,
                       unit_center=units, unit_radius=units)
    elif reg.shape == "ellipse":
        im.mask_ellipse(center=center,
                        radius=[reg.width, reg.height], posangle=pa,
                       unit_center=units, unit_radius=units)
    elif reg.shape == "box":
        im.mask_region(center=center,
                       radius=[reg.width/2.0, reg.height/2.0],
                       posangle=pa, unit_center=units, unit_radius=units)

def _mask_ds9regions(im, regions):
    """Mask regions of a specified MPDAF image, using circle, ellipse
    and/or box region specifications from a ds9 region file.

    Note that only circle, ellipse and box regions are used from the
    list of regions. Other region types are quietly ignored.

    Ds9 region files support various coordinate systems for region
    descriptions. In this function only fk5 coordinates and physical
    (pixel index) coordinates are supported.

    In ds9 a region can be marked to be included or excluded. This
    function goes through the list of regions twice. If any regions
    are marked to be included (the default in ds9), then everything
    outside of those regions is masked. Areas that were already
    masked, and ds9 regions that are marked to be excluded, are also
    masked.

    Parameters
    ----------
    regions : filename or iterable or ds9regions.Ds9Regions
        Either a string that contains the name of a ds9 region file,
        an iterable that returns successive lines of a ds9 file, or
        a Ds9Regions object that contains a parsed list of regions.
    """

    if not isinstance(regions, ds9regions.Ds9Regions):
        try:
            regions = ds9regions.Ds9Regions(regions)
        except ValueError as e:
            raise UserError("Region parsing error: %s" % str(e))
        except IOError as e:
            raise UserError("Region file error, %s" % str(e))

    # Mask all regions that are marked for exclusion.
    for reg in regions:
        if reg.exclude:
            _mask_inside_ds9region(im, reg)

    # Keep a record of the current mask while we determine a
    # separate mask of regions to be kept.
    exclusion_mask = im.mask

    # Accumulate a temporary mask of regions that are to be kept.
    im.mask = np.zeros(im.shape, dtype=bool)
    for reg in regions:
        if not reg.exclude:
            _mask_inside_ds9region(im, reg)

    # Obtain the union of explicitly excluded pixels and
    # pixels outside regions that were explicitly included.
    if im.mask.sum() > 0:
        im.mask = exclusion_mask | ~im.mask
    else:
        im.mask = exclusion_mask

def _plot_fitted_image_results(muse, imfit, muse_im, muse_ft, hst_im, hst_ft,
                               fig=None, display=True, plotfile=None,
                               nowait=False, title=None):
    """Generate a single-page plot of images of the results.

    The plot contains six images. The top left image is the MUSE
    image after it was masked and convolved with the PSF and
    decimation filter of the HST image. The top center image is
    the HST image after it was masked and convolved with the
    best-fit MUSE PSF. The top-right image shows the residual
    image after subtracting these processed MUSE and HST
    images. Below each of these images, the corresponding Fourier
    transform images are shown.

    Parameters
    ----------
    muse : mpdaf.obj.Imagea
       The original MUSE image.
    imfit : FittedImagePhotometry
       The fitted parameter values and other details.
    muse_im : numpy.ndarray
       The fitted version of the MUSE image.
    muse_ft : numpy.ndarray
       The FFT of muse_im.
    hst_im : numpy.ndarray
       The fitted version of the HST image.
    hst_ft : numpy.ndarray
       The FFT of hst_im.
    fig : matplotlib.figure or None
       The figure in which to plot. If this is None, then
       the default figure, plt.gcf(), is substituted.
    display : bool
       If True (the default), display the image.
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
       A specific plot title, or None to request the default title.
       Specify "" if no title is wanted.

    """

    # If no figure has been specified, get the default figure.

    if fig is None:
        fig = plt.gcf()

    # Buffer subsequent graphics until ready to display them.

    plt.ioff()

    # Clear the figure.

    fig.clf()

    # Substitute a default title?

    if title is None:
        title = "MUSE image: %s" % imfit.name

    # Display a plot title?

    if title != "":
        title_y = 0.98 if "\n" in title else 0.95
        fig.suptitle(title, ha="left", x=0.12, y=title_y, fontsize=14)

    # Create a plot grid with 2 rows and 3 columns.

    gs = gridspec.GridSpec(2,3)

    # Determine a suitable range for displaying the the absolute
    # value of the FFT.

    ft_mean = abs(muse_ft).mean()
    ft_std = abs(muse_ft).std()
    ft_vmin = ft_mean - 0.75 * ft_std
    ft_vmax = ft_mean + 4.0 * ft_std

    # Determine a suitable range for displaying the images.

    im_mean = muse_im.mean()
    im_std = muse_im.std()
    im_vmin = im_mean - 2.0 * im_std
    im_vmax = im_mean + 5.0 * im_std

    # Display the MUSE image and its FFT.

    m_im_ax = fig.add_subplot(gs[0,0])
    _plot_2d_array(muse_im, m_im_ax, vmin=im_vmin, vmax=im_vmax,
                  title="MUSE image", ylabel="Pixels",
                  title_fs=13)
    m_ft_ax = fig.add_subplot(gs[1,0])
    _plot_2d_array(abs(muse_ft), m_ft_ax, vmin=ft_vmin, vmax=ft_vmax,
                  title="MUSE FFT", title_fs=13,
                  ylabel="Pixels", xlabel="Pixels")

    # Display the HST image and its FFT.

    h_im_ax = fig.add_subplot(gs[0,1])
    _plot_2d_array(hst_im, h_im_ax, vmin=im_vmin, vmax=im_vmax,
                  title="HST image", title_fs=13)
    h_ft_ax = fig.add_subplot(gs[1,1])
    _plot_2d_array(abs(hst_ft), h_ft_ax, vmin=ft_vmin, vmax=ft_vmax,
                  title="HST FFT", title_fs=13,
                  xlabel="Pixels")

    # Display the difference, MUSE-HST, and its FFT.

    d_im_ax = fig.add_subplot(gs[0,2])
    _plot_2d_array(muse_im - hst_im, d_im_ax,
                  vmin=im_vmin, vmax=im_vmax, title="(MUSE - HST) image",
                  title_fs=13)
    d_ft_ax = fig.add_subplot(gs[1,2])
    _plot_2d_array(abs(muse_ft - hst_ft), d_ft_ax,
                  vmin=ft_vmin, vmax=ft_vmax,
                  title="FFT of (MUSE - HST) image", title_fs=13,
                  xlabel="Pixels")

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

def _save_fitted_images(muse, muse_im, muse_ft, hst, hst_im, hst_ft,
                        crop_indexes, prefix=None):
    """Save the best-fit images from `fit_image_photometry()` to
    FITS files.

    Parameters
    ----------
    muse : mpdaf.obj.Image
       The original MUSE image.
    muse_im : numpy.ndarray
       The fitted version of the MUSE image.
    muse_ft : numpy.ndarray
       The FFT of muse_im.
    hst : mpdaf.obj.Image
       The original HST image.
    hst_im : numpy.ndarray
       The fitted version of the HST image.
    hst_ft : numpy.ndarray
       The FFT of hst_im.
    crop_indexes : [slice, slice]
       The slice indexes that were used to extract subimages
       from the MUSE and HST images.
    prefix : str or None
       The filename prefix to which to append the distinguishing
       suffixes of each image file. If this is None, the original
       filename of the MUSE image is used, after removing the .fits.
    """

    # Get the basename of the file without the FITS suffix.

    if prefix is None:
        prefix = basename(muse.filename).replace(".fits","")

    # Write the images to FITS files.

    _write_image_to_fits(muse, crop_indexes, muse_im,
                         prefix + "_image.fits")
    _write_image_to_fits(hst, crop_indexes, hst_im,
                         prefix + "_hst_model.fits")
    _write_image_to_fits(muse, crop_indexes, muse_im - hst_im,
                         prefix + "_residual.fits")

    # Write the absolute values of the image FFTs to FITS files.

    _write_fft_to_fits(muse, muse_ft, prefix + "_fft.fits")
    _write_fft_to_fits(hst, hst_ft, prefix + "_hst_model_fft.fits")
    _write_fft_to_fits(muse, muse_ft - hst_ft, prefix + "_residual_fft.fits")

def _write_image_to_fits(template, crop_indexes, data, filename):
    """Write an image to a simple FITS file, using
    a specified MPDAF image as a template.

    Parameters
    ----------
    template : `mpdaf.obj.Image`
       The MPDAF image from which the image originated.
    crop_indexes : [slice, slice]
       The slice indexes that were used to extract the
       image from the template before it was processed.
    data : `numpy.ndarray`
       The 2D image array to be written to the file.
    filename : str
       The filename to give the FITS file.
    """

    # Create a primary HDU for the image.

    hdu = fits.PrimaryHDU(data)

    # Get a cropped WCS object to represent the coordinate axes of
    # the image.

    wcs = template.wcs[crop_indexes]

    # Get header keywords for the cropped wcs information.

    wcs_header = wcs.to_header()

    # Add the WCS keywords to the header of the primary HDU.

    hdu.header.extend(wcs_header)

    # Copy selected keywords from the template image.

    for hdr in [template.primary_header, template.data_header]:
        if hdr is not None:
            for key in ['ORIGIN', 'TELESCOP', 'MJD-OBS', 'DATE-OBS', 'OBSERVER',
                        'OBJECT', 'DATE', 'BUNIT', 'FILTER']:
                if key in hdr:
                    hdu.header[key] = (hdr[key], hdr.comments[key])

    # Save the file.

    hdu.writeto(filename, clobber=True)

def _write_fft_to_fits(template, data, filename):
    """Write the absolute values of an FFT of a 2D image to a
    simple FITS file.

    Parameters
    ----------
    template : `mpdaf.obj.Image`
       The MPDAF image from which the FFT was obtained.
    data : `numpy.ndarray`
       The 2D fft array to be written to the file.
    filename : str
       The filename to give the FITS file.
    """

    # Create a primary HDU for the image.

    hdu = fits.PrimaryHDU(np.abs(data))

    # Copy selected keywords from the template image.

    # Copy selected keywords from the template image.

    for hdr in [template.primary_header, template.data_header]:
        if hdr is not None:
            for key in ['ORIGIN', 'TELESCOP', 'MJD-OBS', 'DATE-OBS', 'OBSERVER',
                        'OBJECT', 'DATE', 'BUNIT', 'FILTER']:
                if key in hdr:
                    hdu.header[key] = (hdr[key], hdr.comments[key])

    # Calculate the spatial-frequency increment in cycles per degree.

    im_steps = template.get_step(unit=u.deg)
    ft_steps = 1.0 / im_steps / np.asarray(data.shape)

    # Add header keywords to describe the coordinates.

    hdu.header['CRPIX1'] = data.shape[1] // 2 + 1
    hdu.header['CRPIX2'] = data.shape[0] // 2 + 1
    hdu.header['CRVAL1'] = (0.0, 'The x-axis spatial-frequency origin')
    hdu.header['CRVAL2'] = (0.0, 'The y-axis spatial-frequency origin')
    hdu.header['CD1_1'] = (ft_steps[1], 'cycles/degree')
    hdu.header['CD1_2'] = 0.0
    hdu.header['CD2_1'] = 0.0
    hdu.header['CD2_2'] = (ft_steps[0], 'cycles/degree')

    # Save the file.

    hdu.writeto(filename, clobber=True)

def _generate_fft_images(mfft, hfft):
    """Convert half-plane FFTs of the the best-fit MUSE and HST images
    to full-sized FFTs with zero frequency at the center of the image,
    [ny//2, nx//2].

    Parameters
    ----------
    mfft  : numpy.ndarray(dtype=complex)
       The half-plane FFT of the final MUSE image, as generated by
       np.fft.rfft2().
    hfft  : numpy.ndarray(dtype=complex)
       The half-plane FFT of the final HST image, as generated by
       np.fft.rfft2().

    Returns
    -------
    out : (muse_ft, hst_ft)
       The full complex FFTs of the best-fit MUSE and HST images,
       centered at the center of the image (ie at [ny//2, nx//2]).

    """

    # For display purposes, convert positive half of the conjugate
    # symmetric FFT of the MUSE image to a full FFT and shift
    # frequency 0,0 to pixel index [ny//2, nx//2].

    nx = (mfft.shape[1] - 1) * 2
    ny = mfft.shape[0]
    tmp = np.empty((ny, nx), dtype=mfft.dtype)
    tmp[ny//2:, nx//2:] = mfft[:ny//2, :-1]
    tmp[:ny//2, nx//2:] = mfft[ny//2:, :-1]
    tmp[:ny//2+1, :nx//2] = np.conj(mfft[ny//2::-1, nx//2:0:-1])
    tmp[ny//2+1:,  :nx//2] = np.conj(mfft[ny-1:ny//2:-1, nx//2:0:-1])
    muse_ft = tmp

    # Do the same for the HST FFT.

    nx = (hfft.shape[1] - 1) * 2
    ny = hfft.shape[0]
    tmp = np.empty((ny, nx), dtype=hfft.dtype)
    tmp[ny//2:, nx//2:] = hfft[:ny//2, :-1]
    tmp[:ny//2, nx//2:] = hfft[ny//2:, :-1]
    tmp[:ny//2+1, :nx//2] = np.conj(hfft[ny//2::-1, nx//2:0:-1])
    tmp[ny//2+1:,  :nx//2] = np.conj(hfft[ny-1:ny//2:-1, nx//2:0:-1])
    hst_ft = tmp

    return (muse_ft, hst_ft)

class FitImagePhotometryMP(_FitPhotometryMP):
    """A multiprocessing iterator that creates a pool or worker
    processes to repeatedly call `fit_image_photometry()` for each
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
    kwargs : dict
       An optional dictionary of keyword/value arguments to be passed to
       `fit_image_photometry()`. The default is the empty dictionary, {}.
    nworker : int
       The number of worker processes to use. The default is
       0, which creates multiprocessing.cpu_count() processes.
       Alternatively, if a negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.

    """
    def __init__(self, hst_filename, muse_filenames, kwargs={}, nworker=0):

        # Instantiate the parent _FitPhotometryMP object.

        super(FitImagePhotometryMP, self).__init__(
            hst_filename, muse_filenames, cmd_fn=fit_image_photometry,
            cmd_kwargs=kwargs, nworker=nworker)

    def next(self):
        """Return the results from the next image in the list of input MUSE
        images.

        Returns
        -------
        out : `FittedImagePhotometry`
           The fitting results from the next image.

        """

        return super(FitImagePhotometryMP, self).next()
