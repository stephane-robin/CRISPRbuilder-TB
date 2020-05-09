# CRISPRbuilder-TB

Please refer to the notebook crisprbuilder-tb.ipynb included in the package to see detailed explanations of the package.

## Purpose of this project

Collect and annotate Mycobacterium tuberculosis whole genome sequencing data for CRISPR investigations.

## Requirements

CRISPRbuilder-TB needs the following dependencies to work:

python = "^3.7"
xlrd = "^1.2.0"
openpyxl = "^3.0.3"
xmltodict = "^0.12.0"
biopython = "^1.76"
datetime = "^4.3"
parallel-fastq-dump
balstn+

These different versions are automatically downloaded when installing the CRISPRbuilder-TB package.

## Installation

Install the package by writing in the command prompt: pip3 install CRISPRbuilder-TB.

## How to use this package

The most often used instruction for this package is: python3 CRISPRbuilder-TB --collect {SRA_reference}.

See the attached notebook for a comprehensive explanation.

## History

The current version is 1.0
