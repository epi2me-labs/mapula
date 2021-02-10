"""
Mapping Stats setup script

Copyright (c) 2021 by Oxford Nanopore Technologies Ltd.
"""
import os
from setuptools import setup, find_packages

os.environ["GIT_SSL_NO_VERIFY"] = "true"

__version__ = "1.0.0"

setup(
    name="mapula",
    version=__version__,
    author="epi2melabs",
    setup_requires=["pytest-runner"],
    description="Mapula",
    zip_safe=False,
    install_requires=["pysam", "pandas", "aplanat"],
    packages=find_packages(exclude=("tests",)),
    package_data={
        "mapping_stats": ["mapula/tests/data/*"],
    },
    entry_points={"console_scripts": ["mapula = mapula.main:run_main"]},
)
