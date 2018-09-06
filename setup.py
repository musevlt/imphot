from setuptools import setup, find_packages

setup(
    name='imphot',
    version='0.1',
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
