[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/aminoextract/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/aminoextract)
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/AminoExtract)  

# AminoExtract

AminoExtract is a tool that extracts amino acid sequences from nucleotide sequences based on a GFF input file.

AminoExtract is able to filter the genomic features in the input files to make sure your output makes sense and to write the resulting amino acid sequences to either a single file or to individual files for every feature, depending on your inputs as a user.

## Why this tool?
Because sometimes you just want a dedicated tool to do a mundane task. And sometimes it's just simply necessary to, for example, ensure reproducibility, portability or to facilitate long-term maintainability of larger projects.

Instead of copying this one script used for translating and writing amino acids across all projects, we can now point to AminoExtract for just that.

## Installation requirements

AminoExtract requires python 3.10 or later to work.  

Dependencies such as Pandas, Biopython and python-magic are installed during the installation procedure.

## Installation instructions

Installation through conda and pip will be made available soon.

1. Download the latest version of AminoExtract by cloning this repository and navigate to the newly created direcotry.  
Copy and paste the code-snippet below in order to do so.  
```
git clone https://github.com/RIVM-bioinformatics/AminoExtract.git && cd AminoExtract/
```

2. If necessary, create a conda-environment and install the necessary dependencies.  
Copy and paste the code-snippet below in order to do so.  
```
mamba env create -f env.yml && conda activate AminoExtract
```

3. Now install AminoExtract into the conda environment with the following:  
```
pip install .
```

AminoExtract is now installed!  
You can use AminoExtract from anywhere on your system as long as the conda-environment in which it is installed is active.  
You can test if installation was succesful by typing `AminoExtract -v` which should display the installed version.