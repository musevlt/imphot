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
    scripts=[
        'imphot/scripts/fit_photometry',
        'imphot/scripts/make_wideband_image',
        'imphot/scripts/regrid_hst_to_muse',
    ],
    # entry_points={
    #     'console_scripts': [
    #         'fix-icrs = muse_analysis.scripts.fix_icrs:main'
    #     ],
    # },
)
