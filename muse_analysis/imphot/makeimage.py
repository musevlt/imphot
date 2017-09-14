#!/usr/bin/env python

from __future__ import print_function, division

import math
import astropy.units as u
from scipy import integrate, interpolate
from sys import exc_info
from multiprocessing import (Process, Queue, cpu_count)
from signal import (signal, SIGINT, SIG_IGN)
from traceback import format_exception
import numpy as np
from numpy import ma
from mpdaf.obj import (Cube, Image)

from .core import (UserError, WorkerError)

__all__ = ['bandpass_image']


def bandpass_image(cube, wavelengths, sensitivities,
                   unit_wave=u.angstrom, interpolation="linear",
                   truncation_warning=True, nprocess=0):
    """Given a cube of images versus wavelength and the bandpass
    filter-curve of a wide-band monochromatic instrument, extract
    an image from the cube that has the spectral response of the
    monochromatic instrument.

    For example, this can be used to create a MUSE image that has
    the same spectral characteristics as an HST image. The MUSE
    image can then be compared to the HST image without having to
    worry about any differences caused by different spectral
    sensitivities.

    For each spectral plane of the cube, the filter-curve is
    integrated over the width of that spectral plane to obtain a
    weight, w[n]. The output image is then given by the following
    weighted mean::

       output_image = sum(w[n] * cube_image[n]) / sum(w[n])

    In practice, to accomodate masked pixels, the w[n] array is
    expanded into a cube w[n,y,x], and the weights of individual
    masked pixels in the cube are zeroed before the above equation
    is applied.

    If the wavelength axis of the cube only partly overlaps the
    bandpass of the filter-curve, the filter curve is truncated to
    fit within the bounds of the wavelength axis. A warning is
    printed to stderr if this occurs, because this results in an
    image that lacks flux from some of the wavelengths of the
    requested bandpass.

    This function is equivalent to Cube.bandpass_image() in MPDAF, but
    it is much faster, because of its use of multiple processes.

    Parameters
    ----------
    cube : mpdaf.obj.Cube
      The cube to be processed.
    wavelengths : numpy.ndarray
      An array of the wavelengths of the filter curve,
      listed in ascending order of wavelength. Outside
      the listed wavelengths the filter-curve is assumed
      to be zero.
    sensitivities : numpy.ndarray
      The relative flux sensitivities at the wavelengths
      in the wavelengths array. These sensititivies will be
      normalized, so only their relative values are important.
    unit_wave : `astropy.units.Unit`
      The units used in the array of wavelengths. The default is
      angstroms. To specify pixel units, pass None.
    interpolation : str
      The form of interpolation to use to integrate over the
      filter curve. This should be one of::

        "linear"     : Linear interpolation
        "cubic"      : Cubic spline interpolation (very slow)

      The default is linear interpolation. If the filter curve
      is well sampled and its sampling interval is narrower than
      the wavelength pixels of the cube, then this should be
      sufficient. Alternatively, if the sampling interval is
      significantly wider than the wavelength pixels of the
      cube, then cubic interpolation should be used instead.
      Beware that cubic interpolation is much slower than linear
      interpolation.
    truncation_warning : bool
      By default, a warning is reported if the filter bandpass exceeds
      the wavelength range of the cube. The warning can be disabled
      by setting this argument to False.
    nprocess : int
       The number of worker processes to use. The default is
       0, which creates multiprocessing.cpu_count() processes.
       Alternatively, if a negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.

    Returns
    -------
    out : `~mpdaf.obj.Image`
        An image formed from the filter-weighted mean of spectral
        planes in the cube that overlap the bandpass of the filter
        curve.

    """

    # Where needed, convert the wavelengths and sensitivities
    # sequences into numpy arrays.

    wavelengths = np.asarray(wavelengths, dtype=float)
    sensitivities = np.asarray(sensitivities, dtype=float)

    # The sensitivities and wavelengths arrays must be one
    # dimensional and have the same length.

    if (wavelengths.ndim != 1 or sensitivities.ndim != 1 or
            len(wavelengths) != len(sensitivities)):
        raise UserError('The wavelengths and sensititivies arguments'
                        ' should be 1D arrays of equal length')

    # Convert the array of wavelengths to the wavelength units of the cube.

    if unit_wave != cube.wave.unit:
        if unit_wave is None:
            wavelengths = cube.wave.coord(wavelengths)
        else:
            wavelengths = (wavelengths * unit_wave).to(cube.wave.unit).value

    # Obtain the starting and ending wavelengths of each pixel of the
    # cube.

    indexes = np.arange(cube.shape[0])
    lowave = cube.wave.coord(indexes - 0.5)
    hiwave = cube.wave.coord(indexes + 0.5)
    del indexes

    # Get the indexes of the pixels that contain the starting and
    # ending wavelengths of the filter bandpass.

    start = max(0, int(math.floor(cube.wave.pixel(wavelengths[0]) + 0.5)))
    stop = min(cube.shape[0],
               int(math.floor(cube.wave.pixel(wavelengths[-1]) + 0.5)) + 1)

    # Abort if there is no overlap between the bandpass filter curve
    # and the wavelength coverage of the cube.

    if start >= cube.shape[0] or stop <= 0 or start >= stop:
        raise UserError("The filter curve does not overlap the "
                        "wavelength coverage of the cube.")

    # If the bandpass does not start at the beginning of a wavelength
    # pixel of the cube, modify the effective beginning wavelength for
    # the first overlapping pixel to match it.

    if wavelengths[0] > lowave[start]:
        lowave[start] = wavelengths[0]

    # If the bandpass does not end at the end of a wavelength pixel of
    # the cube, modify the effective ending wavelength of the last
    # overlapping pixel to match it.

    if wavelengths[-1] < hiwave[stop - 1]:
        hiwave[stop - 1] = wavelengths[-1]

    # Obtain an interpolator of the bandpass curve.

    spline = interpolate.interp1d(x=wavelengths, y=sensitivities,
                                  kind=interpolation)

    # If the bandpass curve isn't completely contained within the
    # wavelength range of the cube, tell the user the fraction of the
    # integrated bandpass that lies beyond the bounds of the cube.

    if truncation_warning and \
       (wavelengths[0] < lowave[0] or wavelengths[-1] > hiwave[-1]):

        # Work out the start and stop indexes of the slice needed
        # to truncate the arrays of the bandpass filter curve to
        # the wavelength range of the cube.

        if wavelengths[0] < lowave[0]:
            bpstart = np.searchsorted(wavelengths, lowave[0], 'right')
        else:
            bpstart = 0

        if wavelengths[-1] > hiwave[-1]:
            bpstop = np.searchsorted(wavelengths, hiwave[-1], 'left')
        else:
            bpstop = wavelengths.shape[0]

        # Integrate the overal bandpass filter curve.

        total = integrate.trapz(sensitivities, wavelengths)

        # Also integrate over just the truncated parts of the curve.

        lost = 0.0
        if bpstart > 0:
            s = slice(0, bpstart)
            lost += integrate.trapz(sensitivities[s], wavelengths[s])
        if bpstop < wavelengths.shape[0]:
            s = slice(bpstop, wavelengths.shape[0])
            lost += integrate.trapz(sensitivities[s], wavelengths[s])

        # Compute the fraction of the integrated bandpass response
        # that has been truncated.

        lossage = lost / total

        # Report the loss if it is over 0.5%.

        if lossage > 0.005:
            cube._logger.warning(
                "%.2g%% of the integrated " % (lossage * 100.0) +
                "filter curve is beyond the edges of the cube.")

    # If only one process has been requested, perform the computation
    # here.

    if nprocess == 1:
        sums = _SpectralPlaneSums.combine_spectral_planes(
            cube=cube, start=start, stop=stop, lowave=lowave,
            hiwave=hiwave, spline=spline)

    # Use multiple processes to combine successive ranges of spectral planes?

    else:
        with _WidebandImageMP(cube=cube, start=start, stop=stop,
                              lowave=lowave, hiwave=hiwave,
                              spline=spline, max_planes_per_job=500,
                              nprocess=nprocess) as mp:
            mp.run()
        sums = mp.sums

    # Compute the filter-curve weighted mean data and variance images

    data, var = sums.combined_data_var()

    # Encapsulate the data and variance images in an MPDAF Image
    # object, and return this.

    return Image.new_from_obj(cube, data=data, var=var)


class _SpectralPlaneSums(object):

    '''A class that accumulates the weighted sums of images from multiple
    spectral planes of a MUSE cube, weighted by a filter throughput curve.

    Parameters
    ----------
    cube : `mpdaf.obj.Cube`
       The cube whose images are being summed.
    wsum : `numpy.ndarray`
       Default = None
       An image of the sum of the weights of each pixel. The weight of
       each pixel is the integrated value of the filter curve over the
       wavelength range of its spectral plane, or zero if the pixel is
       masked.

       If this argument is passed None (the default), then arrays of
       zeros are subsituted for the wsum, wdsum and w2vsum arguments.
    wdsum : `numpy.ma.MaskedArray`
       Default = None
       An image of the sum of the flux of each pixel scaled by its weight.
       This is ignored if wsum is None.
    w2vsum : `numpy.ma.MaskedArray` or None
       Default = None
       An image of the sum of the variance of each pixel scaled by its
       weight squared.
       This is ignored if wsum is None.

    Attributes
    ----------
    wsum : `numpy.ndarray`
       An image of the sum of the weights of each pixel. The weight of
       each pixel is the integrated value of the filter curve over the
       wavelength range of its spectral plane, or zero if the pixel is
       masked.
    wdsum : `numpy.ma.MaskedArray`
       An image of the sum of the flux of each pixel scaled by its weight.
    w2vsum : `numpy.ma.MaskedArray` or None
       An image of the sum of the variance of each pixel scaled by its
       weight squared.

    '''

    def __init__(self, cube, wsum=None, wdsum=None, w2vsum=None):

        # Get the shape of the images in the cube.

        image_shape = cube.shape[1:]

        # Create the arrays in which the intermediate image sums will
        # be accumulated. Using masked arrays for all except the sum
        # of weights.

        if wsum is None:
            self.wsum = np.zeros(image_shape)             # Sum[weight]
            self.wdsum = ma.array(np.zeros(image_shape),  # Sum[data * weight]
                                  mask=np.ones(image_shape, dtype=bool))
            self.w2vsum = ma.array(np.zeros(image_shape),  # Sum[var * weight^2]
                                   mask=np.ones(image_shape, dtype=bool))
        else:
            self.wsum = wsum
            self.wdsum = wdsum
            self.w2vsum = w2vsum

    def accumulate_sums(self, sums):
        """Add the image sums from combine_spectral_planes() to the
        sums held by self.

        Parameters
        ----------
        sums : _SpectralPlaneSums
           The sums to be added to those in self.
        """

        # The weight array is not a masked array, and should not
        # contain any invalid values, so we can simply sum successive
        # contributions as a normal numpy array.

        self.wsum += sums.wsum

        # The weighted sums of the data and variances are masked
        # arrays.  Unfortunately, summing two numpy masked arrays
        # masks any elements that are masked in either of the arrays,
        # whereas we want to sum just the unmasked values, and only
        # mask the sums of elements in the two arrays that are both
        # masked. The masked arrays returned by _weighted_image_sums()
        # are returned by np.ma.sum(), which zeros the underlying
        # values of masked elements. This means that we can add the
        # underlying values of masked elements without affecting the
        # final sums, then take the logical AND of the masked values
        # to determine whether to mask the sums.

        self.wdsum.data[:] += sums.wdsum.data
        self.wdsum.mask[:] &= sums.wdsum.mask

        if sums.w2vsum is not None:
            self.w2vsum.data[:] += sums.w2vsum.data
            self.w2vsum.mask[:] &= sums.w2vsum.mask
        else:
            self.w2vsum = None

    def combined_data_var(self):
        '''Return the filter-curve weighted mean data and variance images
        by dividing the weighted sums by the sum of their weights.

        Returns
        -------
        out : (`numpy.ma.MaskedArray`, `numpy.ma.MaskedArray` or None)
           The final image and either the corresponding image of pixel
           variances, or None if the input cube had no variance
           information.
        '''

        data = self.wdsum / self.wsum
        if self.w2vsum is None:
            var = False
        else:
            var = self.w2vsum / self.wsum**2
        return (data, var)

    @classmethod
    def combine_spectral_planes(cls, cube, start, stop, lowave, hiwave,
                                spline):
        """For a specified range of spectral planes of a cube, and a spline
        filter curve, return filter-weighted image sums that can later
        either be combined with other weighted sums, or be divided by
        the sum of weights to yield the final combined data and
        variance images.

        Parameters
        ----------
        cube : mpdaf.obj.Cube
           The cube to be processed.
        start : int
           The index of the first spectral plane to be combined.
        stop : int
           The index of the last spectral plane, plus one.
        spline :  scipy.interpolate.interp1d
           The spline to use to interpolate the filter throughput
           curve, indexed by spectral plane index.

        Returns
        -------
        out : _SpectralPlaneSums
           An object containing the weighted sums.

        """

        # Integrate the bandpass over the range of each spectral pixel
        # to determine the weights of each pixel. For the moment skip
        # the first and last pixels, which need special treatment.
        # Integer pixel indexes refer to the centers of pixels,
        # so for integer pixel index k, we need to integrate from
        # k-0.5 to k+0.5.

        w = np.empty((stop - start))
        for k in range(start, stop):
            w[k - start], err = integrate.quad(spline, lowave[k], hiwave[k])

        # Obtain a sub-cube of the selected spectral planes.

        subcube = cube[start:stop, :, :]

        # Temporarily assume there are no variances.

        var = None
        w2vsum = None

        # Obtain masked arrays of the data and variances for the sums.

        if subcube._mask is ma.nomask:

            # Create a mask tha masks all invalid values in both data and var.

            mask = ~np.isfinite(subcube._data)
            if subcube._var is not None:
                mask |= ~np.isfinite(subcube._var)

            # Create masked arrays that use the above mask.

            data = ma.array(subcube._data, mask=mask)
            if subcube._var is not None:
                var = ma.array(subcube._var, mask=mask)
        else:
            data = subcube.data
            if subcube._var is not None:
                var = subcube.var

        # Create an array of weights for each pixel of the sub-cube, with
        # masked weights zero, and unmasked values equal to the weight of
        # the spectral plane that they belong to.

        wcube = w[:, np.newaxis, np.newaxis] * ~data.mask

        # Calculate the weighted sums for the pixels of the sub-cube.

        wsum = wcube.sum(axis=0)
        wdsum = ma.sum(data * wcube, axis=0)
        if var is not None:
            w2vsum = ma.sum(var * wcube**2, axis=0)

        return cls(cube, wsum, wdsum, w2vsum)


class _WidebandImageWP(Process):
    """A worker process for the WidebandImageMP class.

    Parameters
    ----------
    job_queue : multiprocessing.Queue
       The queue from which to receive job requests from the
       parent process.
    result_queue : multiprocessing.Queue
       The queue to send the results of a job back to the
       parent process.
    cube : mpdaf.obj.Cube
       The MUSE cube to be processed.
    lowave : np.ndarray
       The starting wavelengths of each pixel along the wavelength
       axis of the cube.
    hiwave : np.ndarray
       The ending wavelengths of each pixel along the wavelength
       axis of the cube.
    spline : scipy.interpolate.interp1d
       The spline to use to interpolate the filter throughput
       curve, indexed by spectral plane index.

    """

    def __init__(self, job_queue, result_queue, cube, lowave, hiwave, spline,
                 *args, **kwargs):

        # Instantiate the super-class Process object.

        super(self.__class__, self).__init__(*args, **kwargs)

        self.job_queue = job_queue
        self.result_queue = result_queue
        self.cube = cube
        self.lowave = lowave
        self.hiwave = hiwave
        self.spline = spline

    def run(self):

        # Disable keyboard interrupts, because multiprocessing does
        # something with them that prevents us from catching the
        # KeyboardInterrupt exception, and sometimes they hang the
        # worker processes. We'll let the parent process catch this
        # signal and terminate the worker processes.

        signal(SIGINT, SIG_IGN)

        # Catch all exceptions while running, and forward them
        # to the parent process.

        try:

            # Run the jobs that the parent process sends us.

            while True:

                # Get the next job from the parent process.

                (job, start, stop) = self.job_queue.get()

                # If the parent process has told us to stop, do so.

                if job is None:
                    return

                # Combine the specified range of spectral planes.

                sums = _SpectralPlaneSums.combine_spectral_planes(
                    self.cube, start=start, stop=stop, lowave=self.lowave,
                    hiwave=self.hiwave, spline=self.spline)
                results = [sums, start, stop]
                self.result_queue.put((job, "results", results))

        except Exception as e:
            exc_type, exc_value, exc_traceback = exc_info()
            tb = format_exception(exc_type, exc_value, exc_traceback)
            exc_traceback = None
            self.result_queue.put((0, "exception", str(exc_type), str(e), tb))
            exit(1)


class _WidebandImageMP(object):
    """A multiprocessing object that creates a pool or worker
    processes to repeatedly call combine_spectral_planes() for successive
    ranges of spectral planes, and accumulates the combined image and its
    variances.

    Parameters
    ----------
    cube : mpdaf.obj.Cube
       The cube to be processed.
    start : int
       The first wavelength pixel of the cube that overlaps with the
       bandpass.
    stop : int
       The pixel after the last wavelength pixel of the cube that overlaps
       with the bandpass.
    lowave : np.ndarray
       The starting wavelengths of each pixel along the wavelength
       axis of the cube.
    hiwave : np.ndarray
       The ending wavelengths of each pixel along the wavelength
       axis of the cube.
    spline : scipy.interpolate.interp1d
       The spline to use to interpolate the filter throughput
       curve, indexed by wavelength.
    max_planes_per_job : int
       The maximum number of spectral planes to process at a time in
       each worker process. Lower values reduce the memory usage of
       the script, but at the expense of longer execution
       times. Larger values may result in the script attempting to use
       more memory than is available.
    nprocess : int
       The number of processes to use. The default is 0, which creates
       multiprocessing.cpu_count() processes.  Alternatively, if a
       negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.

    Attributes
    ----------
    sums : _SpectralPlaneSums
       Weighted sums of multiple spectral plane images and their
       variances (if any).

    """

    def __init__(self, cube, start, stop, lowave, hiwave, spline,
                 max_planes_per_job, nprocess=0):

        # Record the cube that is being processed.

        self.cube = cube

        # Record the maximum number of spectral planes per job, the
        # total number of spectral planes, and the number of spectral
        # planes that have been processed so far.

        self.max_planes_per_job = max_planes_per_job

        # The worker processes assume that they are integrating over
        # complete pixels, so set the range of pixels that they cover
        # to all except the first and last pixels, which may only be
        # partly within the bandpass. Note that integer pixel indexes
        # refer to the center of the pixel.

        self.start = start
        self.stop = stop
        self.nplane = self.stop - self.start

        # Keep a record of the total number of jobs, the number of
        # jobs that have been sent to workers, and the number of jobs
        # whose results have been returned by the iterator.

        self.sent_jobs = 0
        self.reported_jobs = 0
        self.running = False

        # The shape of the final image.

        self.image_shape = self.cube.shape[1:]

        # Create the arrays in which the intermediate image sums will
        # be accumulated.

        self.sums = _SpectralPlaneSums(self.cube)

        # If the number of workers is given as zero or a negative
        # number, substitute (number_of_cpus - abs(nprocess)), making
        # sure that this requests at least one worker process.

        if nprocess <= 0:
            nprocess = max(cpu_count() + nprocess, 1)

        # How many spectral planes should be processed in each job?
        # If possible divide the planes between processes such that
        # each process only has to do one job each.

        self.planes_per_job = min(int(math.ceil(float(self.nplane) / nprocess)),
                                  self.max_planes_per_job)

        # How many jobs will be needed?

        self.njob = (self.nplane + self.planes_per_job - 1) // self.planes_per_job

        # If there are no jobs to be done, don't create any
        # multiprocessing resources. The next call to next()
        # will handle this case by returning StopIteration().

        if self.njob < 1:
            return

        # Don't create more worker processes than the number of jobs
        # to be done.

        if nprocess > self.njob:
            nprocess = self.njob

        # Create a queue for sending jobs to the worker processes.

        self.job_queue = Queue()

        # Create a communications queue for receiving results
        # from the worker processes.

        self.result_queue = Queue()

        # Create a list of worker processes.

        self.workers = [_WidebandImageWP(job_queue=self.job_queue,
                                         result_queue=self.result_queue,
                                         cube=cube, lowave=lowave,
                                         hiwave=hiwave, spline=spline)
                        for i in range(nprocess)]

    def __enter__(self):
        """Start the worker processes. This function is called when
        a python context-manager block is entered."""

        # Start the worker processes.

        for worker in self.workers:
            worker.start()

        # Pass one job to each worker, plus an extra job to ensure
        # that there is a pending job waiting when the first one
        # completes.

        for i in range(len(self.workers) + 1):
            self._send_next_job()

        # Set a flag to inform __iter__() and next() that the processes
        # are now running.

        self.running = True
        return self

    def __exit__(self, type, value, traceback):
        """Stop the worker processes. This function is called when
        a python context-manager block exits."""

        # Forcibly terminate any worker processes that haven't already
        # exited.

        for worker in self.workers:
            if worker.is_alive():
                worker.terminate()

        # Get the exit statuses of the worker processes so that these
        # processes don't become zombies. Since join() raise an
        # exception if the worker hasn't been started yet, protect it
        # with a try/except block.

        for worker in self.workers:
            try:
                worker.join()
            except:
                pass

        # Close the communication queues.

        self.job_queue.close()
        self.result_queue.close()

    # Send the next job, if any remain, to the worker processes.

    def _send_next_job(self):

        # If any spectral planes remain to be processed, send the next
        # range of spectral planes to the next free worker process.

        if self.sent_jobs < self.njob:

            # Get the index of the next job.

            job = self.sent_jobs

            # Get the start and stop planes of this job.

            start = self.start + job * self.planes_per_job
            stop = min(start + self.planes_per_job, self.stop)

            # Send the job to the next idle worker.

            self.job_queue.put((job, start, stop))
            self.sent_jobs += 1

    # Combine the exposures.

    def run(self):
        """Use multiple processes to combine multiple exposures of a
        field into a single cube.
        """

        # Is the object ready?

        assert self.running, "Each %s instance must be used within a context manager" % self.__class__.__name__

        # Repeatedly send job requests and receive responses until all
        # all jobs have been finished and reported.

        while self.reported_jobs < self.njob:

            # Get the results of the next job to be completed.

            r = self.result_queue.get()

            # Get the sequential number and status of the job.

            job = r[0]
            status = r[1]

            # Did an exception occur in the worker process?

            if status == "exception":

                # Make sure that any subsequent call to next()
                # returns StopIteration, just in case the caller
                # catches the following exceptions and tries to
                # continue the iteration.

                self.njob = 0

                # Get the details of the exception that occurred
                # in the worker process.

                exc_type = r[2]
                message = r[3]
                traceback = r[4]

                # Re-raise keyboard interrupts, because a
                # keyboard interrupt on the parent only seems to be
                # seen in the child.

                if exc_type == str(KeyboardInterrupt):
                    raise KeyboardInterrupt("Worker process interrupted")

                # Propagate all other unexpected exceptions using our
                # WorkerError that passes on a formatted traceback of the
                # original exception.

                raise WorkerError(message, traceback)

            # No exception occured, so get the latest results.

            (sums, start, stop) = r[2]

            # If any jobs remain to be sent to workers, send a new one
            # to replace the one that just finished.

            self._send_next_job()

            # Add the latest results to the sums of all spectral planes.

            self.sums.accumulate_sums(sums)

            # Record the completion of this job.

            self.reported_jobs += 1

        # Send all workers a stop command.

        for worker in self.workers:
            self.job_queue.put((0, None))

        # Get the exit statuses of the worker processes, so that
        # they don't become zombie processes.

        for worker in self.workers:
            worker.join()
