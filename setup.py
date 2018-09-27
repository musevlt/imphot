import os
from setuptools import setup, find_packages

pkgmeta = {}
pkgmeta_file = os.path.join(os.path.dirname(__file__), 'imphot', 'version.py')
with open(pkgmeta_file) as f:
    code = compile(f.read(), 'version.py', 'exec')
    exec(code, pkgmeta)

setup(
    name='imphot',
    version=pkgmeta['__version__'],
    packages=find_packages(),
    zip_safe=False,
    install_requires=['mpdaf', 'lmfit'],
    entry_points={
        'console_scripts': [
            'fit_photometry = imphot.programs:fit_photometry_main',
            'make_wideband_image = imphot.programs:make_wideband_image_main',
            'regrid_hst_to_muse = imphot.programs:regrid_hst_to_muse_main',
        ],
    },
)
