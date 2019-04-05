from setuptools import setup, find_packages

setup(
    name='imphot',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
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
