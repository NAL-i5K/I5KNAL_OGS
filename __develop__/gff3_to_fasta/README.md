# gff3_to_fasta
The software is used to extract biological sequences (such as spliced transcripts, cds, or peptides) from specific regions of genome based on `GFF3` file.
* Free software: [license](https://github.com/NAL-i5K/I5KNAL_OGS/blob/I5KNAL_OGS/LICENCE.md)

## Features
* **Incoporation of [gff3.py](https://github.com/hotdogee/gff3-py)**: `gff3.py` is contributed by [Han Li](https://github.com/hotdogee) which uses simple data structures to parse a `GFF3` file into a structure composed of simple python `dict` and `list`.
* **Validation**: Validate the [`GFF3` formatting errors](https://github.com/NAL-i5K/I5KNAL_OGS/wiki/QC-phase) utilizing [QC methods](https://github.com/NAL-i5K/I5KNAL_OGS/blob/I5KNAL_OGS/bin/gff-QC.py) contribued by [I5K Workspace@NAL team](https://i5k.nal.usda.gov/). Provide `WARNING` messages for gene models that may have incorrect biological seuqences generated because of GFF3 formatting errors.
* **Easy extraction for biological sequences**: Provide options for extracting six types of biological seuqences.
    - `gene`: gene sequence for each record in the fasta output
        - Gene or pseudogene features need to be included in the gff file
    - `exon`: exon sequence for each record in the fasta output
        - Exon features need to be included in the gff file
    - `g_seq`: genomic sequence (exon and intron regions included) for each record in the fasta output
        - Transcript-level features (such as mRNA, rRNA, pseudogenic transcripts) need to be inlcuded in the gff file
    - `trans`: spliced transcript (only exons included) for each record in the fasta output
        - Exon features are mainly used for splicing. CDS features are used insteand if exon features are absent. If both of cds and exon features are absent, the transcript is not generated and a `WARNING` message is showed with the transcript ID.
    - `cds`: coding sequence (utr exons and introns excluded) for each record in the fasta output
        - CDS features need to be included in the gff file.
    - `pep`: translated peptide sequences (translation based on cds regions) for each record in the fasta output
        - CDS features need to be included in the gff file.
* **Python method for universal translation**: The `translator` method is feasible for 
    - translation from 64 combitions of codons
    - translation from codons with IUB Depiction
    - translation from mRNA (U contained) or CDS (T, instead of U contained)
