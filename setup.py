import sys

from AminoExtract import __prog__, __version__

if sys.version_info.major != 3 or sys.version_info.minor < 10:
    print("Error: you must execute setup.py using Python 3.10 or later")
    sys.exit(1)

from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name=__prog__,
    version=__version__,
    author="Florian Zwagemaker",
    project_urls={"Source Code": "https://github.com/RIVM-bioinformatics/AminoExtract"},
    license="MIT",
    author_email="ids-bioinformatics@rivm.nl",
    long_description=long_description,
    description="AminoExtract is an application to extract aminoacid sequences from a fasta file based on a GFF.",
    long_description_content_type="text/markdown",
    python_requires=">=3.10",
    packages=find_packages(),
    install_requires=["biopython>=1.79", "pandas", "rich==13.*", "python-magic==0.4.*"],
    entry_points={
        "console_scripts": [
            "aminoextract = AminoExtract.__main__:main",
            "AminoExtract = AminoExtract.__main__:main",
        ]
    },
    keywords=[],
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
    ],
)
