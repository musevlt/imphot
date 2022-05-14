import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 5):
    raise Exception('python 3.5 or newer is required')

setup(
    name='imphot',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=find_packages(),
    zip_safe=False,
    include_package_data=True,
    python_requires='>=3.5',
    install_requires=['mpdaf', 'lmfit'],
    extras_require={
        'docs': ['sphinx', 'sphinx_rtd_theme', 'sphinx_automodapi',
                 'numpydoc', 'matplotlib'],
    },
    entry_points={
        'console_scripts': [
            'fit_photometry = imphot.programs:fit_photometry_main',
            'make_wideband_image = imphot.programs:make_wideband_image_main',
            'regrid_hst_to_muse = imphot.programs:regrid_hst_to_muse_main',
        ],
    },
)
