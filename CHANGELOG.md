# Changelog

## [0.3.1](https://github.com/RIVM-bioinformatics/AminoExtract/compare/v0.3.0...v0.3.1) (2023-09-12)


### Bug Fixes

* add the split_attributes keyword to the exposed `get_feature_name_attribute()` function ([f1f6ab7](https://github.com/RIVM-bioinformatics/AminoExtract/commit/f1f6ab704c2304a83c2a4087516563c21be21a0f))

## [0.3.0](https://github.com/RIVM-bioinformatics/AminoExtract/compare/v0.2.1...v0.3.0) (2023-09-11)


### Features

* add `split_attributes` as a keyword to the `GffDataFrame` class and `read_gff()` function to either split or not split the attributes column based on user input ([51ea186](https://github.com/RIVM-bioinformatics/AminoExtract/commit/51ea1868bd389879213bb5ec762cb1b43c4bc56d))
* expose most of the gff handling functions as direct python calls ([fe6a19f](https://github.com/RIVM-bioinformatics/AminoExtract/commit/fe6a19f22c7e250227f48c4bf2532bea870c5407))
* expose the `read_fasta()` function for direct python calls ([b4de7a9](https://github.com/RIVM-bioinformatics/AminoExtract/commit/b4de7a9c7c05fa6fe6c5d262beb736c0cf14282d))


### Bug Fixes

* change file writing to append mode so records won't get overwritten in `writer.py` ([e5f0b0a](https://github.com/RIVM-bioinformatics/AminoExtract/commit/e5f0b0ae524f71869f7fb03173f007fb4d02b9db))
* do not interpret default NA values from strings in pd.read_csv() ([83f821e](https://github.com/RIVM-bioinformatics/AminoExtract/commit/83f821e655507ee83d6b0103a24c87c15b695bb4))
* lift the `_split_attribute_colums()` function out of the GffDataFrame class so it can be used without the GffDataFrame object ([6308ee1](https://github.com/RIVM-bioinformatics/AminoExtract/commit/6308ee132b95eeaa74da8e05cc506889a67d6611))
* solve incorrect file extension when a file is presented with multiple dots ([c349ced](https://github.com/RIVM-bioinformatics/AminoExtract/commit/c349ced06fe4a3a823c4a242771e2791399a876b))


### Dependencies

* add lenient versions to conda environment file and setup.py ([6ccd2c7](https://github.com/RIVM-bioinformatics/AminoExtract/commit/6ccd2c734e872c2aefd526216e62d0286765969b))


### Documentation

* improve the docstrings of several functions in reader.py ([36ab67c](https://github.com/RIVM-bioinformatics/AminoExtract/commit/36ab67c99e4c1df900e26240f356b593664c898d))
* update docstring for `read_gff()` function ([51ea186](https://github.com/RIVM-bioinformatics/AminoExtract/commit/51ea1868bd389879213bb5ec762cb1b43c4bc56d))
* update readme to be more up-to-date regarding some of the newer functionalities ([e305ecb](https://github.com/RIVM-bioinformatics/AminoExtract/commit/e305ecb461bcc2cfe3b500e5f0a3623fdc1c925a))

## [0.2.1](https://github.com/RIVM-bioinformatics/AminoExtract/compare/v0.2.0...v0.2.1) (2022-12-06)


### Bug Fixes

* replace gaps in nucleotide-seq with "N" characters to ensure valid translation when `--keep-gaps` flag is given. (forces ambigious AA call) ([bf4ed5f](https://github.com/RIVM-bioinformatics/AminoExtract/commit/bf4ed5f1492bfc357fe3d64c175c2f7a55e595ee))

## [0.2.0](https://github.com/RIVM-bioinformatics/AminoExtract/compare/v0.1.0...v0.2.0) (2022-12-05)


### Features

* add option to keep gaps in nucleotide sequence during aminoacid translation ([3c2daa8](https://github.com/RIVM-bioinformatics/AminoExtract/commit/3c2daa85db8aff39a6b3bfbcb4fc8d901654a7f7))


### Documentation

* update function docstrings ([fea7643](https://github.com/RIVM-bioinformatics/AminoExtract/commit/fea76438f0e6500926c8db28b01d2c01675219a5))
