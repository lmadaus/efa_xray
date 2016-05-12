from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Setup function
setup(
    name        = 'efa_xray',
    version     = '0.5',
    description = """A simple EnKF implementation for Ensemble Forecast Adjustment
                        (EFA; Madaus and Hakim 2015) using the python-xray
                        library as an ensemble state container """,

    author      = 'Luke Madaus',
    author_email= 'lmadaus@atmos.washington.edu',
    license     = 'MIT',
    classifiers = [
                    'Development Status :: 3 - Alpha',
                    'Topic :: Ensemble Forecasting :: Kalman Filter',
                    'License :: OSI Approved :: MIT License',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3',
                  ],
    packages    = find_packages(exclude['contrib', 'docs', 'tests*']),
    install_requres = ['xarray','numpy','matplotlib','pandas','dask'],
)



                    
