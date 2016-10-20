#!/usr/bin/env python

from __future__ import print_function

import sys
from os.path import basename
import argparse
import astropy.units as u
from astropy.io import fits
import numpy as np
import re
from textwrap import dedent
from mpdaf.obj import (Image, Cube)
from .core import *
from .fitimage import *
from .fitstar import *
from .fitboth import *
from .makeimage import *

__all__ = ['fit_photometry_main', 'regrid_hst_to_muse_main',
           'make_wideband_image_main']

# Create the main function of the fit_photometry program.

def fit_photometry_main(argv):
    """This is the main function of the fit_photometry script.
    It attempts to determine the photometry characteristics of
    a MUSE image via comparisons with an HST image, according to
    arguments passed by the user.

    Parameters
    ----------
    argv : list
       The list of command line arguments. In this list, argv[0] must
       be the name of the calling program. The remaining elements
       should be the arguments. To see what those arguments can be,
       pass argv[1]='--help'.
    """

    # Get the name of the script.

    prog = basename(argv[0])

    # Get a parser for the command-line arguments.

    parser = ImphotArgumentParser(FIT_BOTH, prog)

    # Add optional arguments for selecting how much is printed for each
    # result.

    parser.add_argument('--verbose', action='store_true',
                        help=dedent('''\
                        Report details of each fit, including chi-squared,
                        correlations, etc. Normally only summaries of the
                        fitted parameters are displayed.'''))

    # Add mandatory arguments for the HST image and the MUSE images.

    parser.add_argument('hst_image',
                        help=dedent('''\
                        The name of the FITS file containing the reference
                        HST image. This must have been decimated and
                        resampled onto the same coordinate grid as the MUSE
                        images. The HST bandpass filter used to obtain this
                        image should match the bandpass of the MUSE images as
                        closely as possible. Note that the `regrid_hst_to_muse`
                        script or the `mpdaf.obj.Image.align_with_image()`
                        function can be used create a subimage of an HST
                        image that has the same coordinate grid as a MUSE
                        image.'''))

    parser.add_argument('muse_images', nargs="+",
                        help=dedent('''\
                        One or more FITS files containing MUSE images.  These
                        should all be images of the same field and bandpass
                        as the HST image named by the --hst argument. Note
                        that the `make_wideband_image` script or the
                        `mpdaf.obj.Cube.bandpass_image()` can be used to
                        extract an image of a given HST bandpass from a MUSE
                        cube.'''))


    # Parse the command-line options. Note that this will abort
    # the program if any argument errors are detected.

    options = parser.parse_args(argv[1:])

    # Extract a dictionary of the keyword/value arguments that are
    # specific to fit_image_photometry().

    image_kwargs = extract_function_args(options, fit_image_photometry)

    # If the --star option has been specified, also extract a dictionary of
    # the keyword/value arguments that are specific to fit_image_photometry().

    if options.star is not None:
        star_kwargs = extract_function_args(options, fit_star_photometry)
    else:
        star_kwargs = None

    # True for the first reported line.

    first = True

    # Catch keyboard interrupts and user input errors.

    try:

        # If there is only one MUSE image to be processed, perform the
        # job inline.

        if len(options.muse_images) == 1:

            # If no star fit has been requested, just perform
            # an image fit.

            if options.star is None:
                kwargs = extract_function_args(options, fit_image_photometry)
                imfit = fit_image_photometry(options.hst_image,
                                             options.muse_images[0], **kwargs)
                _report_results(imfit, True, options)

            # If a star fit has been requested, perform both the
            # image and the star fits.

            else:
                kwargs = extract_function_args(options,
                                              fit_image_and_star_photometry)
                (imfit, starfit) = fit_image_and_star_photometry(
                    options.hst_image, options.muse_images[0], **kwargs)
                _report_results(imfit, True, options)
                _report_results(starfit, False, options)

        # When processing more than one MUSE image, use multiple processes.
        # Start with the case where we only want image photometry fits.

        elif options.star is None:

            # Get the subset of arguments for fit_image_photometry().

            kwargs = extract_function_args(options, fit_image_photometry)

            # Use a multiprocessing iterator for the image fitting procedure.

            with FitImagePhotometryMP(options.hst_image,
                                      options.muse_images, kwargs,
                                      nworker=0) as image_iterator:

                # Iterate over the fits of each of the input files.

                for imfit in image_iterator:
                    _report_results(imfit, first, options)
                    first = False

        # Now handle the case where both image and star photometry fits
        # have been requested.

        else:

            # Get the subset of arguments for fit_image_and_star_photometry().

            kwargs = extract_function_args(options,
                                          fit_image_and_star_photometry)

            with FitImageAndStarPhotometryMP(options.hst_image,
                                             options.muse_images, kwargs,
                                             nworker=0) as dual_iterator:

                # Iterate over fits to each of the input files.

                for imfit, starfit in dual_iterator:
                    _report_results(imfit, first, options)
                    first = False
                    _report_results(starfit, first, options)

    # For errors found in the input data, simply report the error
    # then abort.

    except UserError as e:
        print("%s: %s" % (prog, str(e)))
        exit(1)

    # If somebody deliberately kills the program, just abort the
    # program without displaying a traceback.

    except KeyboardInterrupt:
        exit(1)

def _report_results(results, showheader, options):
    if options.verbose:
        print(results)
    print(results.summary(header=showheader))

# Create the main function of the regrid_hst_to_muse script.

def regrid_hst_to_muse_main(argv):
    """This is the main function of the regrid_hst_to_muse script.

    Given a MUSE image or MUSE cube and one or more HST images with
    30mas pixels, it resamples the HST images onto the pixel grid of
    the MUSE image, and scales their flux units from electrons/s to
    the flux units of the MUSE image (usually 1e-20
    erg/cm^2/s/Angstrom).

    Parameters
    ----------
    argv : list
       The list of command line arguments. In this list, argv[0] must
       be the name of the calling program. The remaining elements
       should be the arguments. To see what those arguments can be,
       pass argv[1]='--help'.

    """

    # Get the name of the script.

    prog = basename(argv[0])

    # Catch keyboard interrupts and user-input errors.

    try:

        # Get a parser for the command-line arguments.

        parser = argparse.ArgumentParser(prog)

        # Specify optional the arguments.

        parser.add_argument('--quiet', action='store_true',
                            help='''
                            Suppress the messages that report each step as it
                            is performed.''')

        parser.add_argument('--field', nargs='?', default="",
                            metavar="name",
                            help='''
                            When this option used, it specifies a field name
                            to use in the output filename instead of using
                            the input filename. The output files are called
                            hst_<filter>_for_<field>.fits, where
                            <filter> is the name of the filter taken
                            from the header of the HST FITS file, and <field>
                            is either the basename of the MUSE input file
                            (minus its .fits extension), or the value of this
                            optional parameter.''')

        # Add mandatory arguments for the MUSE image and the HST images.

        parser.add_argument('muse_image',
                            help='''
                            The filename of a template MUSE image or MUSE
                            cube in FITS format.
                            ''')

        parser.add_argument('hst_images', nargs="+",
                            help='''
                            The filenames of one or more HST images
                            with 30mas pixels.
                            ''')

        # Parse the command-line options. Note that this will abort
        # the program if any argument errors are detected.

        options = parser.parse_args(argv[1:])

        # Read the MUSE image or MUSE cube.

        try:
            if not options.quiet:
                print("Reading MUSE image: %s" % options.muse_image)
            try:
                muse = Image(options.muse_image)
                if muse.ndim != 2:
                    raise ValueError("Invalid image dimensions")
            except:
                muse = Cube(options.muse_image)
                if muse.ndim != 3:
                    raise ValueError("Invalid cube dimensions")

        except Exception as e:
            print("Error reading: %s (%s)" % (options.muse_image, e),
                  file=sys.stderr)
            exit(1)

        # Process the HST images sequentially.

        for hst_image in options.hst_images:

            # Read the HST FITS file.

            try:
                if not options.quiet:
                    print("Reading HST image: %s" % hst_image)

                # If there is no DATA or SCI extension, MPDAF will
                # refuse to read the file unless it is told which
                # extension to read. Tell it to read the primary data
                # array in this case.

                hdulist = fits.open(hst_image)
                ext = None if ('DATA' in hdulist or 'SCI' in hdulist) else 0
                del hdulist

                # Read the HST image.

                hst = Image(hst_image, ext=ext)

            except Exception as e:
                print("Error reading: %s (%s)" % (hst_image, e),
                      file=sys.stderr)
                exit(1)

            # Get the name of the HST filter.

            filter_name = HstFilterInfo(hst).filter_name.upper()

            # Resample the HST image onto the coordinate grid of the MUSE image.
            if not options.quiet:
                print("Resampling the HST image onto the MUSE pixel grid.")
            regrid_hst_like_muse(hst, muse, inplace=True)

            # Rescale the HST image to have the same units as the MUSE image.

            if not options.quiet:
                print("Changing the flux units of the HST image to match the MUSE image")
            rescale_hst_like_muse(hst, muse, inplace=True)

            # Compose the name of the output file.

            if options.field != "":
                field = options.field
            else:
                field = basename(muse.filename).replace(".fits","")
            output = "hst_" + filter_name + "_for_" + field + ".fits"

            # Write the output file.

            if not options.quiet:
                print("Writing the output file: %s" % output)
            hst.write(output, savemask="dq")

    # For errors found in the input data, simply report the error
    # then abort.

    except UserError as e:
        print("%s: %s" % (prog, e))
        exit(1)

    except KeyboardInterrupt:
        exit(1)


# Create the main function of the make_wideband_image script.

def make_wideband_image_main(argv):
    """This is the main function of the make_wideband_image script.

    Given a MUSE cube and the wavelength sensitivity curve of a
    monochromatic camera, extract an image from the cube that has the
    same wavelength sensitivity curve as the camera.

    Parameters
    ----------
    argv : list
       The list of command line arguments. In this list, argv[0] must
       be the name of the calling program. The remaining elements
       should be the arguments. To see what those arguments can be,
       pass argv[1]='--help'.

    """

    # Get the name of the script.

    prog = basename(argv[0])

    # Catch keyboard interrupts and user-input errors.

    try:

        # Get a parser for the command-line arguments.

        parser = argparse.ArgumentParser(prog)

        # Specify optional the arguments.

        parser.add_argument('--prefix', nargs='?', default="",
                            metavar="output-prefix",
                            help='''
                            When this option used, it specifies a prefix to
                            use when constructing the names of the output
                            FITS files.  The output files are then named,
                            "<prefix>_<input_name>", where <prefix> is the
                            specified string, and <input_name> is the
                            base-name of the input cube, but with the first
                            instance of either, "datacube" or "cube",
                            replaced with "image". When this option is
                            omitted, the prefix defaults to the base-name of
                            the filter file.''')

        parser.add_argument('--quiet', action='store_true',
                            help='''
                            Suppress the messages that report each step as it
                            is performed.''')

        parser.add_argument('--wave_units', nargs='?', default=u.angstrom,
                            metavar="wavelength-units",
                            type=u.Unit,
                            help='''DEFAULT=%(default)s.
                            The units of the wavelengths in the file of
                            sensitivity versus wavelength. If this
                            option is not specified, angstroms are
                            assumed.''')

        parser.add_argument('--nprocess', nargs='?', default=0,
                            metavar="number-of-processes",
                            type=int,
                            help='''DEFAULT=%(default)s.
                            The number of worker processes to use. The
                            default is 0, which creates one process
                            per CPU.  Alternatively, if a negative
                            number, -n, is specified, then max(ncpu -
                            n,1) processes are created, where ncpu is
                            the number of available CPUs. To prevent
                            the use of any worker processes, pass 1 to
                            this argument, and the computation will be
                            performed entirely within the current
                            process.''')

        # Add mandatory arguments for the MUSE image and the HST images.

        parser.add_argument('filter_curve',
                            help='''
                            The filename of the wavelength sensitivity
                            curve of the camera. The file should be a
                            text file with two columns of numbers. The
                            first column are the filter- curve
                            wavelengths, and the second are the
                            corresponding sensitivities, which should be
                            unit-less. The sensitivities will be
                            normalized by the area under the
                            sensitivity curve, so only their relative
                            values are important.
                            ''')

        parser.add_argument('muse_cubes', nargs="+",
                            help='''The filenames of a MUSE
                            cubes to be processed.''')

        # Parse the command-line options. Note that this will abort
        # the program if any argument errors are detected.

        options = parser.parse_args(argv[1:])

        # Read the filter response curve.

        try:
            if not options.quiet:
                print("Reading filter: %s" % options.filter_curve)
                curve = np.loadtxt(options.filter_curve, usecols=(0,1))
        except Exception as e:
            print("Error reading: %s (%s)" % (options.filter_curve, e),
                  file=sys.stderr)
            exit(1)

        # Get the prefix to use for the output files.

        if options.prefix != "":
            prefix = options.prefix
        else:
            prefix = basename(options.filter_curve).split('.')[0]

        # Process the cubes sequentially.

        for muse_cube in options.muse_cubes:

            # Read the cube.

            try:
                if not options.quiet:
                    print("Reading cube: %s" % muse_cube)
                cube = Cube(muse_cube)
            except Exception as e:
                print("Error reading: %s (%s)" % (muse_cube, e),
                      file=sys.stderr)
                exit(1)

            # Compute the output image.

            if not options.quiet:
                print("Computing the output image.")
            image = bandpass_image(cube, curve[:,0], curve[:,1],
                                   unit_wave=options.wave_units,
                                   nprocess=options.nprocess)

            # Get the input filename.

            name = basename(muse_cube)

            # Perform a case-insensitive search for the first occurrence of
            # the word "datacube", or "cube"

            m = re.search(r'(data)?cube', name, re.IGNORECASE)

            # If a match was found, replace it with the word "image",
            # using the case of the first letter of the match.

            if m:
                r = "IMAGE" if m.group(0)[0].isupper() else "image"
                name = name.replace(m.group(0), r, 1)

            # Create the filename of the output file.

            output = prefix + '_' + name

            # Save the resulting image.

            if not options.quiet:
                print("Writing image to: %s" % output)
            image.write(output, savemask="dq")

    # For errors found in the input data, simply report the error
    # then abort.

    except UserError as e:
        print("%s: %s" % (prog, e))
        exit(1)

    except KeyboardInterrupt:
        exit(1)
