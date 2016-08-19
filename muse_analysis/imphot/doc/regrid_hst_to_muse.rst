.. _regrid_hst_to_muse:

Resampling HST images onto the pixel grid of a MUSE image
=========================================================

Both of the photometry fitting algorithms depend on having an HST
image that is sampled on the same pixel grid as the MUSE image that is
being characterized. A suitable HST image can be obtained by applying
a script called :ref:`regrid_hst_to_muse<regrid_hst_to_muse>` to a higher
resolution HST image. The :ref:`regrid_hst_to_muse<regrid_hst_to_muse>` script
resamples the high-resolution HST image onto the same pixel grid as
the MUSE image, and scales its pixel fluxes from HST electrons s\
:sup:`-1` to the flux units of the MUSE image (usually 1e-20 erg cm\
:sup:`-2` s\ :sup:`-1` Angstrom\ :sup:`-1`).

The arguments of the function can be obtained by running it with the
`-h` option::


  % regrid_hst_to_muse -h
  usage: regrid_hst_to_muse [-h] [--quiet] [--field [name]]
                            muse_image hst_images [hst_images ...]

  positional arguments:
    muse_image      The filename of a template MUSE image in FITS format.
    hst_images      The filenames of one or more HST images with 30mas pixels.

  optional arguments:
    -h, --help      show this help message and exit
    --quiet         Suppress the messages that report each step as it is
                    performed.
    --field [name]  When this option used, it specifies a field name to use in
                    the output filename instead of using the input filename. The
                    output files are called hst_<filter>_for_<field>.fits, where
                    <filter> is the name of the filter taken from the header of
                    the HST FITS file, and <field> is either the basename of the
                    MUSE input file (minus its .fits extension), or the value of
                    this optional parameter.
  %

The :ref:`regrid_hst_to_muse<regrid_hst_to_muse>` script gives each resampled
output file a name, ``hst_<filter_name>_for_<muse_filename>``, where
``<filter_name>`` is the name of the HST filter found in the FILTER
keyword of the FITS header, and ``<muse_filename>`` is the name of the
input MUSE FITS file, including its ``.fits`` suffix.

In the following example, the script is used to resample HST images
taken through 2 different filters, called F606W and F775W::

  % regrid_hst_to_muse IMA-R-MUSE-2014-09-24T05\:49\:36.365.fits hlsp_xdf_hst_acswfc-30mas_hudf_*.fits
  Reading MUSE image: IMA-R-MUSE-2014-09-24T05:49:36.365.fits
  Reading HST image: hlsp_xdf_hst_acswfc-30mas_hudf_f606w_v1_sci.fits
  Resampling the HST image onto the MUSE pixel grid.
  Changing the flux units of the HST image to match the MUSE image
  Writing the output file: hst_F606W_for_IMA-R-MUSE-2014-09-24T05:49:36.365.fits
  Reading HST image: hlsp_xdf_hst_acswfc-30mas_hudf_f775w_v1_sci.fits
  Resampling the HST image onto the MUSE pixel grid.
  Changing the flux units of the HST image to match the MUSE image
  Writing the output file: hst_F775W_for_IMA-R-MUSE-2014-09-24T05:49:36.365.fits
  %

By default the output filenames are based on the HST filter name and
the name of the input MUSE file. If the MUSE file is just one exposure
of multiple observations of the same field, then the regridded HST
image will be usable for any of them. In this case it may be
less confusing to indicate a name for the field, such that this gets
used in the output filename, rather than the name of one
exposure. This can be done using the optional ``--field`` argument, as
follows::

  % regrid_hst_to_muse --field UDF01 IMA-R-MUSE-2014-09-24T05\:49\:36.365.fits hlsp_xdf_hst_acswfc-30mas_hudf_*.fits
  Reading MUSE image: IMA-R-MUSE-2014-09-24T05:49:36.365.fits
  Reading HST image: hlsp_xdf_hst_acswfc-30mas_hudf_f606w_v1_sci.fits
  WARNING: MpdafUnitsWarning: Error parsing the BUNIT: 'ELECTRONS/S' did not parse as unit: At col 0, ELECTRONS is not a valid unit. Did you mean electron? [mpdaf.obj.data]
  Resampling the HST image onto the MUSE pixel grid.
  Changing the flux units of the HST image to match the MUSE image
  Writing the output file: hst_F606W_for_UDF01.fits
  Reading HST image: hlsp_xdf_hst_acswfc-30mas_hudf_f775w_v1_sci.fits
  Resampling the HST image onto the MUSE pixel grid.
  Changing the flux units of the HST image to match the MUSE image
  Writing the output file: hst_F775W_for_UDF01.fits
  %

Note that the output filenames in this case, with the addition of the
option ``--field UDF01``, changed to ``hst_F606W_for_UDF01.fits`` and
``hst_F775W_for_UDF01.fits``.

Resampling from within a python script
--------------------------------------

The resampling operation can alternatively be performed from within an
arbitrary python script using the `imphot.regrid_hst_like_muse()`
function, followed by `imphot.rescale_hst_like_muse()`. For example::

  from mpdaf.obj import Image
  import imphot

  hst = Image(hst_fits_filename)
  muse = Image(muse_fits_filename)
  imphot.regrid_hst_like_muse(hst, muse, inplace=True)
  imphot.rescale_hst_like_muse(hst, muse, inplace=True)


