import os

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

VERSION = "0.4.2"

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
    version=VERSION,
    packages=find_packages(),
    author="Wiep van der Toorn",
    author_email="w.vandertoorn@fu-berlin.de",
    include_package_data=True,
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
        "tqdm",
        "numpy==1.24.4",
        "pandas",
        "joblib",
    ],
    extras_require={
        "live-demux": [
            "attrs",
            "minknow-api==5.7.2",
        ],
    },
)
