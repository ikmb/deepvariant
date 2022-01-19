# Usage information

## Running the pipeline

This pipeline requires nextflow and singularity, so make sure these modules are loaded.

To run the pipeline, a typical command will looks as follows:

`nextflow run ikmb/deepvariant --samples Samples.csv --genome GRCh38`

These options are further explained in the following:

### Options

#### `--samples`
This option expects a sample sheet in CSV format that specifies information on the subject, sample, details of the library and the location of the paired-end files.

A script is included that takes a folder full of fastQ files and creates a compliant sample sheet. This script requires ruby and can be downloaded [here](https://github.com/ikmb/deepvariant/blob/master/bin/samplesheet_from_folder.rb)

* Illumina short reads

`ruby samplesheet_from_folder.rb --folder /path/to/reads > Samples.csv`

* Pacbio CCS/HiFi reads

`ruby samplesheet_from_folder.rb --folder /path/to/reads --pacbio > Samples.csv`

Please note that the script makes some assumptions about the naming structure of the FastQ files. For Illumina reads, the script tries to group reads by library name into patient/sample sets. For Pacbio, this is not possible - so you will have to edit the sample sheet to properly set up the grouping of multiple movies into one sample and patient. 

The columns indivID and sampleID can otherwise be changed by you to better suit your needs. The indivID determines to output folder name (i.e. a patient), the sampleID should be designate a particular biological sample. 
### `--pacbio`

Specifies that the input data are pacbio CCS reads. This triggers an alternative analysis workflow, including SV calling.

#### `--genome`
The name of the genome assembly version to use. Allowed options are:

* GRCh37 (1000Genomes reference with decoys)
* GRCh38 (Current human genome WITHOUT alt contigs, as recommended by Heng Li)
* hg38 (Current human genome with all alt contigs, as distributed by the BROAD)

For Illumina input, please note that all assemblies come with pre-defined calling regions, as provided by the BROAD institute. These were manually curated to exclude regions that cannot be reliably called with short read data.
