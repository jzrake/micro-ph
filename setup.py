#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

long_description = """
EOS stuff
"""

setup(
    name="eospy",
    version="0.4.0",
    author="Jonathan Zrake",
    author_email="zrake@nyu.edu",
    packages=["eospy"],
    url="http://github.com/jzrake/micro-ph",
    license="MIT",
    description="Compute microphysics terms in astrophysical EOS.",
    long_description=long_description,
    ext_modules = [Extension('eospy._eospy', ['src/micro-ph.c',
                                              'src/solvers.c',
                                              'eospy/_eospy.c'])],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
