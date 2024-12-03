import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

import re

VERSIONFILE = "warpdemux/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

ext_modules = [
    Extension(
        name=str("warpdemux.segmentation._c_segmentation"),
        sources=[str(os.path.join("warpdemux", "segmentation", "_c_segmentation.pyx"))],
        include_dirs=[np.get_include()],
        language="c++",
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

setup(
    name="WarpDemuX",
    version=verstr,
    packages=find_packages(exclude=["warpdemux/adapted*"]),
    author="Wiep van der Toorn",
    author_email="w.vandertoorn@fu-berlin.de",
    include_package_data=True,
    package_data={
        "warpdemux.config": ["config_files/*.toml"],
        "warpdemux.models": ["model_files/*.toml", "model_files/*.joblib"],
    },
    ext_modules=cythonize(ext_modules, language_level="3"),
    cmdclass={"build_ext": build_ext},
    entry_points={"console_scripts": ["warpdemux = warpdemux.main:main"]},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Typing :: Typed",
        "License :: CC BY-NC 4.0",
    ],
    install_requires=[
        "pod5",
        "vbz-h5py-plugin",
        "dtaidistance",
        "scikit-learn==1.3.1",
        "scipy",
        "cython==0.29.36",
        "toml",
        "torch==2.4.1",
        "tqdm",
        "numpy==1.24.4",
        "pandas",
        "joblib",
        "attrs",
        #include "bottleneck"
        "bottleneck",
    ],
    extras_require={
        "live-demux": [
            "minknow-api==5.7.2",
        ],
    },
)
