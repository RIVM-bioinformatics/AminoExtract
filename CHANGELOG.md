# Changelog

## [0.4.0](https://github.com/RIVM-bioinformatics/AminoExtract/compare/v0.3.1...v0.4.0) (2025-09-17)


### Features

* add a write out to gff method to the GFFDataframe class ([d2fffb2](https://github.com/RIVM-bioinformatics/AminoExtract/commit/d2fffb22a402a44e9ce46c8876fc32785632393d))
* added spliced gene support ([4b597bb](https://github.com/RIVM-bioinformatics/AminoExtract/commit/4b597bbddab4e2236df63c8e12a6c66d3797fc5e))
* AmpliGone.reader - added support for space seperated GFF file instead of only ; separeted, formatted with black ([eac5985](https://github.com/RIVM-bioinformatics/AminoExtract/commit/eac5985e2d296f8a9ecb60a9afaf15c80671f9a2))


### Bug Fixes

* "seqid" column is also fine for splicing info instead of just "ID" column ([d567c54](https://github.com/RIVM-bioinformatics/AminoExtract/commit/d567c54e69f5cc9fd13171cca52e0314f03f0b96))
* add configless gff reader access ([e81bc29](https://github.com/RIVM-bioinformatics/AminoExtract/commit/e81bc29573946d294ab41e3e1d4f66d5b3442b81))
* add header to gff write out ([977a8e7](https://github.com/RIVM-bioinformatics/AminoExtract/commit/977a8e79f9b2bdbc031a2aab640117a75c6b7ff6))
* add Name column to features even if not in gff ([b003a66](https://github.com/RIVM-bioinformatics/AminoExtract/commit/b003a662d842461726553cd84b8c7c6c04aab465))
* added better Name attribute normalization ([56f5638](https://github.com/RIVM-bioinformatics/AminoExtract/commit/56f56381dd8485cf9bcc6be77902ddb5f1e27104))
* also added splicing information for samples without splicing. Fixed 1-based inclusiveness. ([17a6ed9](https://github.com/RIVM-bioinformatics/AminoExtract/commit/17a6ed91e54f9cea580ce1a2ee0e5844151f5544))
* also fix when using get_feat_name_attr ([037fe59](https://github.com/RIVM-bioinformatics/AminoExtract/commit/037fe5977d33c78f01d92e6f6c888fe3da29a3b5))
* AminoExtract.logging - fixed typing ([4d3f446](https://github.com/RIVM-bioinformatics/AminoExtract/commit/4d3f4465ae3f25633a8aa9eee3e0ab4cbc5a547c))
* AminoExtract.reader - added _normalize_attributes step to change gene_name into Name ([7b4c1fe](https://github.com/RIVM-bioinformatics/AminoExtract/commit/7b4c1fe55f7a4a98b4447499b017078da6469e6b))
* complex get_feature_name_attributes with many differnt feature_types ([17f0645](https://github.com/RIVM-bioinformatics/AminoExtract/commit/17f06453e410ace7fe6e8ee6c924fd3240352f79))
* env.yml - intel channel is no longer accesible, throws an error ([eac5985](https://github.com/RIVM-bioinformatics/AminoExtract/commit/eac5985e2d296f8a9ecb60a9afaf15c80671f9a2))
* fixed all pylint issues ([3aa58f2](https://github.com/RIVM-bioinformatics/AminoExtract/commit/3aa58f24a945dce926295ac6262b9a64d8d5c01e))
* fixed issues and security hotspots identified by sonarqube ([457239d](https://github.com/RIVM-bioinformatics/AminoExtract/commit/457239d5b94aa5a46de3f7ab50ec974df5dc47de))
* fixed string ending with ; issue ([9e130fc](https://github.com/RIVM-bioinformatics/AminoExtract/commit/9e130fc0438deaed13006b3858c0f3e29b19f1e1))
* include optional seqid instead of ID ([72765d8](https://github.com/RIVM-bioinformatics/AminoExtract/commit/72765d85e351c7af01cbaf12d7cb197101a9c8fc))
* remove colnames from gff ([a818684](https://github.com/RIVM-bioinformatics/AminoExtract/commit/a818684be081a15575dda4e5560e5eab65e15df1))
* removed logging because there is no log ([bbf44b4](https://github.com/RIVM-bioinformatics/AminoExtract/commit/bbf44b4df3b6d0f01edaa1f88da443e7138b027a))
* split different CDS into new lines in complex GFF files ([4797735](https://github.com/RIVM-bioinformatics/AminoExtract/commit/4797735c091c49584dd48543cedbea90532a863d))
* update GffDataFrame dtype definition to support int-based seqids ([eac306f](https://github.com/RIVM-bioinformatics/AminoExtract/commit/eac306f26d828ebb5aa0b92020eb940623086f22))
* using intel as a channel causes the installation to fail due to a HTTP 403 (Forbidden) error. As it is not required for the depencies it can safely be removed. ([9056c1f](https://github.com/RIVM-bioinformatics/AminoExtract/commit/9056c1fa54b5701237ff52212d655d9b51af29a1))

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
