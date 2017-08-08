# solarbextrapolation
[![Build Status](https://travis-ci.org/sunpy/solarbextrapolation.svg?branch=master)](https://travis-ci.org/sunpy/solarbextrapolation)
[![Documentation Status](http://readthedocs.org/projects/solarbextrapolation/badge/?version=latest)](http://docs.sunpy.org/projects/solarbextrapolation/en/latest/?badge=latest)


solarbextrapolation is a library for extrapolating 3D magnetic fields from line-of-sight magnetograms. 

## Installation
First, install the following dependencies using your preferred Python package manager. We recommend [the Anaconda distribution](https://www.continuum.io/downloads).

* numPy
* matplotlib
* sciPy
* astropy
* sunPy

Optionally, if you'd like to do advanced visualizations of the extrapolation field, you'll need the mayavi package.

solarbextrapolation is **not** yet available on PyPI or conda-forge. To install the package,

```shell
$ git clone https://github.com/sunpy/solarbextrapolation.git
$ cd solarbextrapolation
$ python setup.py install
```

and to run the tests,

```shell
$ python setup.py test
```

