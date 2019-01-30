from sys import exc_info
from multiprocessing import (Process, Queue, cpu_count)
from signal import (signal, SIGINT, SIG_IGN)
from traceback import format_exception

from mpdaf.obj import Image

from .core import (UserError, WorkerError)


class _FitPhotometryWP(Process):
    """A worker process for the FitPhotometryMP class.

    Parameters
    ----------
    job_queue : multiprocessing.Queue
       The queue from which to receive job requests from the
       parent process.
    result_queue : multiprocessing.Queue
       The queue to send the results of a job back to the
       parent process.
    hst_filename : str
       The name of the HST image FITS file.
    cmd_fn : function
       The fitting function to be called. This must be compatible
       with the call signature, cmd_fn(hst, muse, **kwargs),
       where hst and muse are mpdaf.obj.Image objects containing
       HST and MUSE images.
    cmd_kwargs : dict
       A dictionary of the keyword-argument values to pass to
       fit_image_photometry().

    """

    def __init__(self, job_queue, result_queue, hst_filename, cmd_fn,
                 cmd_kwargs={}, *args, **kwargs):

        # Instantiate the super-class Process object.

        super(_FitPhotometryWP, self).__init__(*args, **kwargs)

        self.job_queue = job_queue
        self.result_queue = result_queue
        self.hst_filename = hst_filename
        self.cmd_fn = cmd_fn
        self.cmd_kwargs = cmd_kwargs

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

            # Attempt to read the HST image.

            try:
                hst = Image(self.hst_filename)

                # If an exception occurs, return a formatted traceback of
                # the exception.

            except Exception:
                raise UserError("Error reading HST file: %s" % self.hst_filename)

            # Run the jobs that the parent process sends us.

            while True:

                # Get the next job from the parent process.

                (job, muse_filename) = self.job_queue.get()

                # If the parent process has told us to stop, do so.

                if muse_filename is None:
                    return

                # Attempt to open the MUSE file.

                try:
                    muse = Image(muse_filename)
                except Exception:
                    raise UserError("Error reading MUSE file: %s" % muse_filename)

                # Perform the photometry fitting operation.

                results = self.cmd_fn(hst, muse, **self.cmd_kwargs)
                self.result_queue.put((job, "results", results))

        except Exception as e:
            exc_type, exc_value, exc_traceback = exc_info()
            tb = format_exception(exc_type, exc_value, exc_traceback)
            exc_traceback = None
            self.result_queue.put((0, "exception", str(exc_type), str(e), tb))
            exit(1)


class _FitPhotometryMP(object):
    """A multiprocessing iterator that creates a pool or worker
    processes to repeatedly call a specified function for each
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
    cmd_fn : function
       The fitting function to be called. This must be compatible
       with the call signature, cmd_fn(hst, muse, **kwargs),
       where hst and muse are mpdaf.obj.Image objects containing
       HST and MUSE images.
    cmd_kwargs : dict
       An optional dictionary of keyword/value arguments to be passed to
       cmd_fn(). The default is the empty dictionary, {}.
    nworker : int
       The number of worker processes to use. The default is
       0, which creates multiprocessing.cpu_count() processes.
       Alternatively, if a negative number, -n, is specified, then
       max(multiprocessing.cpu_count()-n,1) processes are created.

    """

    def __init__(self, hst_filename, muse_filenames, cmd_fn, cmd_kwargs={},
                 nworker=0):

        # Keep a record of the names of the MUSE files to be processed.

        self.muse_filenames = muse_filenames

        # Keep a record of the total number of jobs, the number of
        # jobs that have been sent to workers, and the number of jobs
        # whose results have been returned by the iterator.

        self.njob = len(muse_filenames)
        self.sent_jobs = 0
        self.reported_jobs = 0
        self.running = False

        # If there are no jobs to be done, don't create any
        # multiprocessing resources. The next call to next()
        # will handle this case by returning StopIteration().

        if self.njob < 1:
            return

        # Create a list for the results of the jobs.

        self.job_results = [None for i in range(self.njob)]

        # If the number of workers is given as zero or negative
        # number, substitute (number_of_cpus - abs(nworker)), making
        # sure that this requests at least one worker process.

        if nworker <= 0:
            nworker = max(cpu_count() + nworker, 1)

        # Don't create more worker processes than the number of jobs
        # to be done.

        if nworker > self.njob:
            nworker = self.njob

        # Create a queue for sending jobs to the worker processes.

        self.job_queue = Queue()

        # Create a communications queue for receiving results
        # from the worker processes.

        self.result_queue = Queue()

        # Create a list of worker processes.

        self.workers = [_FitPhotometryWP(job_queue=self.job_queue,
                                         result_queue=self.result_queue,
                                         hst_filename=hst_filename,
                                         cmd_fn=cmd_fn,
                                         cmd_kwargs=cmd_kwargs)
                        for i in range(nworker)]

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
        """Start the worker processes. This function is called when
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

        # If any files remain to be processed, send the next to
        # the next free worker process.

        if self.sent_jobs < self.njob:
            job = self.sent_jobs
            self.job_queue.put((job, self.muse_filenames[job]))
            self.sent_jobs += 1

    # Throw an exception if __enter__() hasn't been called.

    def _assert_running(self):
        assert self.running, "Each %s instance must be used within a context manager" % self.__class__.__name__

    # Mark this object as an interator.

    def __iter__(self):

        # Is the iterator ready?

        self._assert_running()

        # Return the iterator.

        return self

    # Return the next value.

    def next(self):
        """Return the results from the next image in the list of input MUSE
        images.

        Returns
        -------
        out : `FittedPhotometry`
           The fitted results from the next image.

        """

        # Is the iterator ready?

        self._assert_running()

        # Stop if all jobs have been finished and reported.

        if self.reported_jobs >= self.njob:

            # Send all workers a stop command.

            for worker in self.workers:
                self.job_queue.put((0, None))

            # Get the exit statuses of the worker processes, so that
            # they don't become zombie processes.

            for worker in self.workers:
                worker.join()

            # Tell the caller that the iterator has finished.

            raise StopIteration()

        # Receive results from worker processes until the results of
        # the next job to be returned by the iterator has been
        # received. Note that self.reported_jobs is the index of the
        # job that we want next.

        while self.job_results[self.reported_jobs] is None:

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

                # Preserve UserError exceptions, so that the caller
                # can respond to errors caused by incorrect user input.

                if exc_type == str(UserError):
                    raise UserError(message)

                # Also re-raise keyboard interrupts, because a
                # keyboard interrupt on the parent only seems to be
                # seen in the child.

                if exc_type == str(KeyboardInterrupt):
                    raise KeyboardInterrupt("Worker process interrupted")

                # Propagate all other unexpected exceptions using our
                # WorkerError that passes on a formatted traceback of the
                # original exception.

                raise WorkerError(message, traceback)

            # No exception occured, so record the results of the fit.

            self.job_results[job] = r[2]

            # If any jobs remain to be sent to workers, send a new one
            # to replace the one that just finished.

            self._send_next_job()

        # The results of the next job of the iterator has been
        # received, so return it.

        result = self.job_results[self.reported_jobs]
        self.reported_jobs += 1
        return result
