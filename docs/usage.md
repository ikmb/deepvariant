# Usage information

## Running the pipeline

This pipeline requires nextflow and singularity, so make sure these modules are loaded.

To run the pipeline, a typical command will looks as follows:

`nextflow run ikmb/deepvariant --samples Samples.csv --genome GRCh38_p14`

These options are further explained in the following:

### Analytical options

#### `--tools`

The pipeline supports different tool chains, depending on input. These can be specified as comma-separated list, for example:

Short reads: `--tools Deepvariant,Manta,Vep` will perform read alignments, SNP/Indel as well as SV calling and finally run the variant effect predictor on the SNP/INDEL set(s).

Long reads: `--tools Deepvariant,PBSV` will perform read alignments and SNP/Indel as well as SV calling.

Providing no tools will only perform alignment and deduplication of alignments. 

All options:

* Deepvariant - perform variant calling (long- and short reads)
* PBSV - perform SV calling (long reads)
* Manta - perform SV calling (short reads)
* VEP - perform variant effect prediction (long- and short reads, requires Deepvariant)
* Intersect - create a set of FastQ files that only include those reads overlapping --intervals (or the built-in calling regions)

### Sequencing setup(s)

WGS data is often produced in multiple sequencing runs or across multiple lanes. This pipeline will automatically deal with such situations by using the sampleID as a grouping variable. Please do not
concatenate your reads beforehand ; just ensure that data belonging to the same biological sample are grouped under the same id. 

### Options

#### `--samples`
This option expects a sample sheet in CSV format that specifies information on the subject, sample, details of the library and the location of the paired-end files.

A script is included that takes a folder full of fastQ files and creates a compliant sample sheet. This script requires ruby and can be downloaded [here](https://github.com/ikmb/deepvariant/blob/master/bin/samplesheet_from_folder.rb)

* Illumina short reads

`ruby samplesheet_from_folder.rb --folder /path/to/reads > Samples.csv`

```
IndivID;SampleID;libraryID;rgID;rgPU;R1;R2
NA24143_I33978;Sample_NA24143_I33978;NA24143_I33978-L2;HHNVKDRXX.2.NA24143_I33978-L2;HHNVKDRXX.2.32142;/path/to/library_R1_001.fastq.gz;/path/to/library_R2_001.fastq.gz
```

* Pacbio CCS/HiFi reads

`ruby samplesheet_from_folder.rb --folder /path/to/reads --pacbio > Samples.csv`

``` 
IndivID;SampleID;R1
Indiv_NA24143_I33978;Sample_NA24143_I33978;/path/to/movie.fastq.gz
```

Please note that the script makes some assumptions about the naming structure of the FastQ files. For Illumina reads, the script tries to group reads by library name into patient/sample sets. For Pacbio, this is not possible - so you will have to edit the sample sheet to properly set up the grouping of multiple movies into one sample and patient. 

The columns indivID and sampleID can otherwise be changed by you to better suit your needs. The indivID determines to output folder name (i.e. a patient), the sampleID should be designate a particular biological sample. 
#### `--pacbio`

Specifies that the input data are pacbio CCS reads. This triggers an alternative analysis workflow, including SV calling.

#### `--genome`
The name of the genome assembly version to use. Allowed options are:

* GRCh37 (1000Genomes reference with decoys)
* GRCh38 (DEFUNCT: human genome WITHOUT alt contigs, as recommended by Heng Li)
* GRCh38_p14 (Current human reference genome recommended for variant analysis)
* hg38 (the Wittig reference, p13 with only the primary chromosomes and no unplaced contigs etc)
* CHM13v2 (the telomere-to-telomere reference with fully assembled Y chromosome)

For Illumina input, please note that all assemblies come with pre-defined calling regions, as provided by the BROAD institute. These were manually curated to exclude regions that cannot be reliably called with short read data. The only exception is the CHM13v2 assembly, for which no such regions are currently available - so run times might be a bit increased. 

#### `--intervals`

A BED file with calling regions; otherwise the whole-genome calling region set from BROAD will be used. 

#### `--phase`

Enable phasing for short-read data (phasing is performed by default for Pacbio data)

