import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def circular_moffat_profile(rsq, fwhm, beta):

    """Return the value of a circularly symmetric Moffat profile at a
    specified value of the radius-squared from the center of the
    profile.

    Parameters
    ----------
    rsq  :  float
       The radius squared from the center of the moffat.
    fwhm   :  float
       The full-width at half-maximum of the Moffat.
    beta   : float
       The exponent of the moffat.

    Returns
    -------
    out : numpy.ndarray
       The value of the Moffat profile at each radius-squared in rsq.

    """

    # A Moffat function is defined as follows:
    #
    #   y(r,a,beta) = 1 / (1 + r**2 / a**2)**beta
    #
    # Calculate a**2 from the FWHM and beta.

    asq = fwhm**2 / 4.0 / (2.0**(1.0 / beta) - 1.0)

    # Compute the Moffat function at each value of rsq.

    return 1.0 / (1.0 + rsq / asq)**beta

def circular_gaussian_profile(rsq, fwhm):

    """Return the value of a circularly symmetric Gaussian profile at a
    specified value of the radius-squared from the center of the
    profile.  The profile is normalized to be 1.0 at rsq=0.0.

    Parameters
    ----------
    rsq  :  float
       The radius squared from the center of the gaussian.
    fwhm   :  float
       The full-width at half-maximum of the Gaussian.

    Returns
    -------
    out : numpy.ndarray
       The value of the Gaussian profile at each radius-squared in rsq.

    """

    # A Gaussian function is defined as follows:
    #
    #   y(x,y) = exp(-0.5 * r**2 / sigma**2)
    #
    # In terms of the FWHM, sigma**2 = fwhm**2  / (8*ln(2)
    #
    # Calculate arg = -0.5/sigma**2

    arg = -4.0 * np.log(2.0)/ fwhm**2

    # Compute the Gaussian function at each value of rsq.

    return np.exp(arg * rsq)

def main():

    # If true, display the plot interactively. If false write it to a PDF file.

    interactive = False

    # Set the number of pixels along the X and Y axes of the image.
    # It's best to make this a power of 2.

    n = 1024

    # Set the FWHM of the Moffat and Gaussian PSFs to be drawn.

    fwhm = 1.0

    # Set the beta value of the Moffat function.

    beta = 2.0

    # Set the width of the image, in the same units as the FWHM of the
    # PSFs.

    width = 50.0

    # Set the maximum radial distance of the profile to be plotted.

    xmax = 2.0    # Arcsec

    # Set the maximum radial frequency of the FFT to be plotted.

    fmax = 1.5    # Cycles per arcsec

    # Calculate the interval between pixels in the image domain.

    dx = width / n

    # Calculate the interval between pixels in the Fourier domain.

    df = 1.0 / width

    # Calculate the coordinate of the center of each pixel along
    # the X and Y axes.

    x = np.fft.fftfreq(n,df)

    # Calculate the frequency coordinate of the center of each Fourier
    # transform pixel along the X and Y axes.

    f = np.fft.fftfreq(n,dx)

    # Calculate the radius squared of each pixel in the image,
    # relative to the origin and wrap-around rules that are assumed by
    # the FFT algorithm.

    rsq = np.fft.fftfreq(n, df)**2 + \
          np.fft.fftfreq(n, df)[np.newaxis,:].T**2

    # Compute the pixels of the Gaussian and Moffat images.

    g_im = circular_gaussian_profile(rsq, fwhm)
    m_im = circular_moffat_profile(rsq, fwhm, beta)

    # Compute the FFTs of the above images.

    g_ft = np.fft.fft2(g_im)
    m_ft = np.fft.fft2(m_im)

    # Create the figure that will host the plots.

    if interactive:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=(7.5,3.5), dpi=600, tight_layout=True)

    # Configure a grid of 2 plots, one above the other.

    gs = gridspec.GridSpec(1, 2)

    # In the first plot, draw a 1D cut through the positive part of
    # the Gaussian and Moffat profiles.

    ax = fig.add_subplot(gs[0,0])
    ax.set_autoscale_on(False)
    ax.set_xlim(0.0, xmax)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel("Radial distance (arcsec)")
    ax.set_ylabel("Fraction of peak")
    xslice = slice(0,int((xmax+dx) / dx)+1)
    gim_trace, = ax.plot(x[xslice], g_im[0, xslice],ls=':')
    mim_trace, = ax.plot(x[xslice], m_im[0, xslice])
    ax.legend([gim_trace, mim_trace], ["Gaussian PSF", "Moffat PSF"],
              loc='upper right', fontsize=10, labelspacing=0.3, frameon=False)
    ax.text(0.6*xmax, 0.5, "Image Plane", ha="center")

    # In the second plot, draw a 1D cut through the positive part of
    # the Fourier transforms of the Gaussian and Moffat profiles.

    ax = fig.add_subplot(gs[0,1])
    ax.set_autoscale_on(False)
    ax.set_xlim(0.0, fmax)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel("Radial frequency (cycles/arcsec)")
    ax.set_ylabel("Fraction of peak")
    fslice = slice(0,int((fmax+df) / df)+1)
    gim_trace, = ax.plot(f[fslice], np.real(g_ft[0,fslice]/g_ft[0,0]),ls=':')
    mim_trace, = ax.plot(f[fslice], np.real(m_ft[0,fslice]/m_ft[0,0]))
    ax.legend([gim_trace, mim_trace], ["Gaussian PSF", "Moffat PSF"],
              loc='upper right', fontsize=10, labelspacing=0.3, frameon=False)
    ax.text(0.6*fmax, 0.5, "Fourier Plane", ha="center")
    if interactive:
        plt.show()
    else:
        fig.savefig("psf_fft.pdf", landscape=True)

main()
