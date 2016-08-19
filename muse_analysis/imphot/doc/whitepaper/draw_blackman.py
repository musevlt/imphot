#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def circular_blackman_profile(f, fc):

    """Return the value of a circularly symmetric Blackman window at a
    specified value of the radius-squared from the center of the
    profile.

    Parameters
    ----------
    f    :  float
       The radial frequency from the center of the profile.
    fc   :  float
       The cutoff frequency of the blackman window (ie. where
       it falls to zero).

    Returns
    -------
    out : numpy.ndarray
       The value of the Blackman profile at each radius-squared in rsq.

    """

    # Relative to its center, the Blackman window function is defined
    # as follows:
    #
    #   if abs(f) <= fc:
    #     y(f,fc) = 0.42 + 0.5 * cos(pi*f/fc) + 0.08 * cos(2*pi*f/fc)
    #   elif abs(f) > fc:
    #     y(f,fc) = 0.0

    theta = np.pi * f / fc
    return np.where(np.abs(theta) <= np.pi,
                    0.42 + 0.5 * np.cos(theta) + 0.08 * np.cos(2*theta), 0.0)

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

    n = 2048

    # Set the cutoff frequency of the blackman window to half the sampling
    # frequency of MUSE images, in cycles per arcsecond.

    fc = 0.5 / 0.2

    # Set the width of the Fourier transform image, in cycles/arcsecond.

    freq_width = fc * 50.0

    # Set the maximum radial distance of the profile to be plotted.

    xmax = 1.5      # Arcsec

    # Set the maximum radial frequency of the FFT to be plotted.

    fmax = 1.5*fc   # Cycles per arcsec

    # Calculate the interval between pixels in the Fourier domain.

    df = freq_width / n

    # Calculate the interval between pixels in the image domain.

    dx = 1 / freq_width

    # Calculate the coordinate of the center of each pixel along
    # the X and Y axes.

    x = np.fft.fftfreq(n,df)

    # Calculate the frequency coordinate of the center of each Fourier
    # transform pixel along the X and Y axes.

    f = np.fft.fftfreq(n,dx)

    # Calculate the radius squared of each pixel in the Fourier
    # transform image, relative to the origin and wrap-around rules
    # that are assumed by the FFT algorithm.

    fsq = f**2 + f[np.newaxis,:].T**2

    # Calculate the radius squared of each pixel in the image plane,
    # relative to the origin and wrap-around rules that are assumed by
    # the FFT algorithm.

    rsq = x**2 + x[np.newaxis,:].T**2

    # Compute the pixels of the Blackman window in the Fourier plane.

    b_ft = circular_blackman_profile(np.sqrt(fsq), fc)

    # Compute the inverse FFT of the above images and normalize it.

    b_im = np.real(np.fft.ifft2(b_ft))
    b_im /= b_im[0,0]

    # Calculate the FWHM of the blackman window.

    #fwhm = 0.810957290099188 * fc

    # Find the FWHM of the inverse Fourier transform of the Blackman window.
    # (This was found empirically)

    psf_fwhm = 1.21638214428 / fc

    # Calculate the FWHM of a Gaussian that when used as a window in the
    # fourier plane, produces a Gaussian PSF of the above FWHM.

    gwin_fwhm = 4.0 * np.log(2.0) / (np.pi * psf_fwhm)

    # Create a 2D image of a circular gaussian of the above FWHM in the
    # Fourier plane, and also the inverse FFT of this image.

    g_ft = circular_gaussian_profile(fsq, gwin_fwhm)
    g_im = np.real(np.fft.ifft2(g_ft))
    g_im /= g_im[0,0]

    # Create the figure that will host the plots.

    if interactive:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=(7.5,3.5), dpi=600)

    # Configure a grid of 2x2 plots.

    gs = gridspec.GridSpec(2, 2)

    # Reduce the gap between the axes of the upper and lower plots.

    gs.update(hspace=0.05, wspace=0.25, bottom=0.15)

    # Draw a 1D cut through the positive part of the Gaussian and
    # Blackman profiles, plotting twice; the second time magnified
    # to show the tails of the functions close to zero.

    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax1.set_autoscale_on(False)
    ax2.set_autoscale_on(False)
    ax1.set_xlim(0.0, xmax)
    ax1.set_xticklabels([])
    ax2.set_xlim(0.0, xmax)
    ax1.set_ylim(-0.1, 1.1)
    ax2.set_ylim(-0.0025, 0.0055)
    set_yticks(ax2, 0.001, 0.002)
    ax2.set_xlabel("Radial distance (arcsec)")

    # Get the slice prescription for extracting a radial cut along the
    # X-axis.

    sx = slice(0,int((xmax+dx) / dx)+1)

    # Get the radii of the slice elements.

    x_slice = x[sx]

    # Plot a slice along the positive X-axis of the Gaussian image.

    g_slice = g_im[0, sx]
    gim1_trace, = ax1.plot(x_slice, g_slice,ls=':')
    gim2_trace, = ax2.plot(x_slice, g_slice,ls=':')

    # Plot a slice along the positive X-axis of the Blackman image.

    b_slice = b_im[0, sx]
    bim1_trace, = ax1.plot(x_slice, b_slice)
    bim2_trace, = ax2.plot(x_slice, b_slice)

    # Plot a legend that identifies the two traces.

    ax1.legend([gim1_trace, bim1_trace], ["Gaussian IFT", "Blackman IFT"],
              loc='upper right', fontsize=10, labelspacing=0.3, frameon=False)
    ax1.text(0.65*xmax, 0.4, "Image Plane", ha="center")

    # Draw a 1D cut through the positive part of the Gaussian and
    # Blackman FFT profiles, plotting twice; the second time magnified
    # to show the tails of the functions close to zero.

    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,1])
    ax1.set_autoscale_on(False)
    ax2.set_autoscale_on(False)
    ax1.set_xlim(0.0, fmax)
    ax1.set_xticklabels([])
    ax2.set_xlim(0.0, fmax)
    ax1.set_ylim(-0.1, 1.1)
    ax2.set_ylim(-0.0025, 0.0055)
    set_yticks(ax2, 0.001, 0.002)
    ax2.set_xlabel("Radial frequency (cycles/arcsec)")

    # Get the slice prescription for extracting a radial cut along the
    # X-axis.

    sx = slice(0,int((fmax+df) / df) + 1)

    # Get the radii of the slice elements.

    f_slice = f[sx]

    # Plot a slice along the positive X-axis of the Gaussian FFT.

    g_slice = g_ft[0, sx]
    gim1_trace, = ax1.plot(f_slice, g_slice,ls=':')
    gim2_trace, = ax2.plot(f_slice, g_slice,ls=':')

    # Plot a slice along the positive X-axis of the Blackman FFT.

    b_slice = b_ft[0, sx]
    bim1_trace, = ax1.plot(f_slice, b_slice)
    bim2_trace, = ax2.plot(f_slice, b_slice)

    # Plot a legend that identifies the two traces.

    ax1.legend([gim1_trace, bim1_trace], ["Gaussian window", "Blackman window"],
              loc='upper right', fontsize=10, labelspacing=0.3, frameon=False)
    ax1.text(0.65*fmax, 0.4, "Fourier Plane", ha="center")
    ax1.plot([fc,fc], [-0.1,0.1])
    ax2.plot([fc,fc], [-0.0025,-0.0018])
    ax1.text(fc, 0.12, r"$f_c$", ha="center", va="bottom", color="red",
             fontsize=14)
    ax2.text(fc, -0.0016, r"$f_c$", ha="center", va="bottom", color="red",
             fontsize=14)

    # Plot a common Y-axis label.

    fig.text(0.02, 0.5, "Fraction of peak", rotation="vertical", ha="left",
             va="center")

    if interactive:
        plt.show()
    else:
        fig.savefig("blackman_ift.pdf", landscape=True)

def set_xticks(ax, minor, major):
    """Change the major and/or minor ticks of an X-axis.

    Parameters
    ----------
    ax    : Axes
       This axes of the plot.
    minor : float or None
       The desired interval between minor ticks, or None to leave
       the minor ticks unchanged.
    major : float or None
       The desired interval between major ticks, or None to leave
       the major ticks unchanged.

    """

    # Get the limits of the axis.

    vmin, vmax = ax.get_xlim()

    # Define the major ticks.

    if major is not None:
        ax.set_xticks(np.arange(np.ceil(vmin/major)*major,
                                np.ceil(vmax/major)*major,major), minor=False)

    # Define the minor ticks.

    if minor is not None:
        ax.set_xticks(np.arange(np.ceil(vmin/minor)*minor,
                                np.ceil(vmax/minor)*minor,minor), minor=True)

def set_yticks(ax, minor, major):
    """Change the major and/or minor ticks of an Y-axis.

    Parameters
    ----------
    ax    : Axes
       This axes of the plot.
    minor : float or None
       The desired interval between minor ticks, or None to leave
       the minor ticks unchanged.
    major : float or None
       The desired interval between major ticks, or None to leave
       the major ticks unchanged.

    """

    # Get the limits of the axis.

    vmin, vmax = ax.get_ylim()

    # Define the major ticks.

    if major is not None:
        ax.set_yticks(np.arange(np.ceil(vmin/major)*major,
                                np.ceil(vmax/major)*major,major), minor=False)

    # Define the minor ticks.

    if minor is not None:
        ax.set_yticks(np.arange(
            np.ceil(vmin/minor), np.ceil(vmax/minor)) * minor, minor=True)

main()
