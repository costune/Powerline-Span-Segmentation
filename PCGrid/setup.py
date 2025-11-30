import os
import sys
import subprocess
import platform
from pathlib import Path

from setuptools import setup, Extension

def get_pybind_extension():
    import pybind11
    
    source_files = [
        "ext.cpp",
        "src/Las2PowerLine.cpp",
        "src/PCGrid.cpp",
        "src/PointCloutFeature.cpp",
    ]
    
    include_dirs = [
        pybind11.get_include(),
        "src",
        "src/nanoflann",
        "third_party/libLAS/include",
        "/usr/include/eigen3",
        "/usr/local/include",
    ]
    
    # library_dirs = [
    #     "/usr/local/lib",
    #     "third_party/libLAS/build/bin/Release"
    # ]
    
    extra_compile_args = [
        "-std=c++11",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Wno-unused-parameter",
        "-Wno-unused-variable",
    ]

    # extra_link_args = [
    #     "-Wl,-rpath,/home/jyxc/projects/PowerLineLidar_tmp/PCGrid/cmake_build/third_party/libLAS/bin/Release"
    # ]

    ext_modules = [
        Extension(
            "pcgrid._C",
            sources=source_files,
            include_dirs=include_dirs,
            # library_dirs=library_dirs,
            # libraries=libraries,
            extra_compile_args=extra_compile_args,
            # extra_link_args=extra_link_args,
            language="c++",
        ),
    ]
    
    return ext_modules



try:
    ext_modules = get_pybind_extension()
    cmdclass = {}
except ImportError:
    print("Warning: pybind11 not found. Please install it: pip install pybind11")
    ext_modules = []
    cmdclass = {}



setup(
    name="pcgrid",
    version="0.0.0",
    packages=["pcgrid"],
    ext_modules=ext_modules,
    cmdclass=cmdclass,
    install_requires=[
        "numpy>=1.18.0",
    ]
)
