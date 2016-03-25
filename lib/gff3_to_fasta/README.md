=============
gff3_to_fasta
=============

The software is used to extract biological sequences (such as spliced transcripts, cds, or peptides) from specific regions of genome based on |GFF3|_ file.
* Free software: [license](https://github.com/NAL-i5K/I5KNAL_OGS/blob/I5KNAL_OGS/LICENCE.md)

Features
--------

* **Incoporation of [gff3.py](https://github.com/hotdogee/gff3-py)**: `gff3.py` is contributed by [Han Li](https://github.com/hotdogee) which uses simple data structures to parse a |GFF3|_ file into a structure composed of simple python `dict` and `list`.
* **Validation**: Validate the |GFF3|_ [formatting errors](https://github.com/NAL-i5K/I5KNAL_OGS/wiki/QC-phase) utilizing [QC methods](https://github.com/NAL-i5K/I5KNAL_OGS/blob/I5KNAL_OGS/bin/gff-QC.py) contribued by [I5K Workspace@NAL team](https://i5k.nal.usda.gov/). Provide `WARNING` messages for gene models that may have incorrect biological seuqences generated because of |GFF3|_ formatting errors.
* **Easy extraction for biological sequences**: Provide options for extracting six types of biological seuqences.
    - **`gene`**: Gene sequence for each record in the fasta output. Gene or pseudogene features need to be included in the gff file
    - **`exon`**: Exon sequence for each record in the fasta output. Exon features need to be included in the gff file
    - **`pre_trans`**: Genomic region of a transcript model, namely premature transcript (exon and intron regions included), for each record in the fasta output. Transcript-level features (such as mRNA, rRNA, pseudogenic transcripts) need to be inlcuded in the gff file.
    - **`trans`**: Spliced transcript (only exons included) for each record in the fasta output. Exon features are mainly used for splicing. CDS features are used instead if exon features are absent. If both of cds and exon features are absent, the transcript is not generated and a `WARNING` message is showed with the transcript ID.
    - **`cds`**: Coding sequence (utr exons and introns excluded) for each record in the fasta output. CDS features need to be included in the gff file.
    - **`pep`**: Translated peptide sequences (translation based on cds regions) for each record in the fasta output. CDS features need to be included in the gff file.
* **`translator` method for universal translation**: The `translator` method is feasible for 
    - translation from 64 combitions of codons
    - translation from codons with IUB Depiction
    - translation from mRNA (U contained) or CDS (T, instead of U contained)

Quick Start
-----------

Sequences extraction based on genome annotation file in |GFF3|_ format might generate incorrect sequences if the |GFF3|_ file contain certain errors. This program can automatically validate |GFF3|_ format and list the detected errors for users along with the extracted seuqences by simple command as follows.

`python gff3_to_fasta.py -g sample_input/annotations3.gff -f sample_input/sample.fa -st cds -d simple -o sample_output`

If you would like to ignore QC step for checking errors in gff file, you can use the bellow command.

`python gff3_to_fasta.py -g sample_input/annotations3.gff -f sample_input/sample.fa -st cds -d simple -o sample_output -noQC`

If you would like to incorporte a specfic method in `gff3_to_fata` into other python script, bollow is a suggested script.

.. code:: python

import gff3_to_fasta

seq = 'GTGGCTCGTTTGATTGAACAAATATGTACTAACCCAGTTGGATTATCTGGATCTGGATTTTTTCTGGTGACAAAGAATTTTCTACTTCAGATGGCAGGAACGATAGTTACATTTGAACTGATGCTGTTTCAATTTGCCCCAGTAAATGCACAGCAAAAACCCATGAAGTCATATGACTGTATTTAA'

gff3_to_fasta.translator(seq)


.. |GFF3| replace:: ``GFF3``
.. |dict| replace:: ``dict``
.. |list| replace:: ``list``
.. |FASTA| replace:: ``FASTA``

.. _GFF3: http://www.sequenceontology.org/gff3.shtml
.. _dict: https://docs.python.org/2/tutorial/datastructures.html#dictionaries
.. _list: https://docs.python.org/2/tutorial/datastructures.html#more-on-lists
.. _FASTA: http://en.wikipedia.org/wiki/FASTA_format

