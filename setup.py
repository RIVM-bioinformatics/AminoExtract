import sys

from AminoExtract import __prog__, __version__

if sys.version_info.major != 3 or sys.version_info.minor < 7:
    print("Error: you must execute setup.py using Python 3.7 or later")
    sys.exit(1)

from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name=__prog__,
    version=__version__,
    author="Florian Zwagemaker",
    project_urls={"Source Code": "https://github.com/florianzwagemaker/AminoExtract"},
    license="MIT",
    author_email="florian.zwagemaker@rivm.nl",
    long_description=long_description,
    description="AminoExtract is a basic program to extra aminoacid sequences from a fasta file based on a GFF.",
    long_description_content_type="text/markdown",
    python_requires=">=3.7",
    packages=find_packages(),
    install_requires=["biopython", "gffpandas"],
    entry_points={
        "console_scripts": [
            "aminoextract = AminoExtract:main",
            "AminoExtract = AminoExtract:main",
        ]
    },
    keywords=[],
    zip_safe=False,
)
