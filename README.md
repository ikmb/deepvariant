![](images/ikmb_bfx_logo.png)

# IKMB DeepVariant Pipeline

This pipeline performs an end-to-end variant calling, starting from raw fastQ files to a multi-sample VCF. It has been pre-configured for the CCGA MedCluster. 

## Running the pipeline

This pipeline requires nextflow and singularity, so make sure these modules are loaded. 

To run the pipeline, a typical command will looks as follows:

`nextflow run ikmb/deepvariant --samples Samples.csv --genome GRCh38` 

These options are further explained in the following:

### Options

#### `--samples` 
This option expects a sample sheet in CSV format that specifies information on the subject, sample, details of the library and the location of the paired-end files. 

A script is included that takes a folder full of fastQ files and creates a compliant sample sheet. This script requires ruby and can be downloaded [here](https://github.com/ikmb/deepvariant/bin/samplesheet_from_folder.rb)

`ruby samplesheet_from_folder.rb --folder /path/to/reads > Samples.csv`

#### `--genome`
The name of the genome assembly version to use. Allowed options are:

* GRCh37 (1000Genomes reference with decoys)
* GRCh38 (Current human genome WITHOUT alt contigs, as recommended by Heng Li)
* hg38 (Current human genome with all alt contigs, as distributed by the BROAD)

Please note that all assemblies come with pre-defined calling regions, as provided by the BROAD institute. These manually were curated to exclude regions that cannot be reliably called with short read data. 


