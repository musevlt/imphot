from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution('imphot').version
except DistributionNotFound:
    # package is not installed
    __version__ = None
