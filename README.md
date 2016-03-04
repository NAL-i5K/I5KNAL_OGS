# I5KNAL_OGS
This project is to develop python tools for generating official gene set (OGS) by integrating manually curated and predicted gene annotations (GFF3 format). There are two phases involved: (1) QC phase and (2) Merge phase. A prototype of the whole pipeline has been done by I5K Workspace@NAL team. However, the source codes of the prototype program is not release for public, because it incorporated several components written by programming languages other than Python. Therefore, this project will re-implement those non-python components, and expects to deliver a complete python package for OGS generation. If you have urgent needs for OGS generation, you can send queries to I5K [at] ars.usda.gov. The i5k team can help to host your data, and apply OGS generation pipeline on your data for you.

## \_\_develop\_\_/
Tools under development.
* example_file/
    - Example files for testing
* function4gff/
    - Functions for dealing with gff3
* gff3_to_fasta/
    - Extract specific sequeces from genome sequences accroding to gff file.
* inter_model/
    - QC functions for processing multiple features between models (inter-model) in GFF3 file.
* intra_model/
    - QC functions for processing multiple features within a model (intra-model) in GFF3 file.
* single_feature/
    - QC functions for processing every single feature in GFF3 file.
* template/
    - Template script for development

## bin/
General script for running through different phases of the OGS pipeline.
* gff-QC.py
    - Detection of GFF format errors (Error_code: Error_tag)
        + Esf0001: pseudogene or not?
        + Esf0002: Negative/Zero start/end coordinate
        + Ema0005: unusual child features in the type of pseudogene found
        + Emr0001: Duplicate transcripts found

## lib/
Completed tools would be shown as under a specific directory. Tools under development would be shown as a Symbolic link.
* gff3_modified/
    - Basic data structure used for nesting the information of genome annotations in GFF3 format.
