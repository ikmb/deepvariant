![](images/ikmb_bfx_logo.png)

# IKMB DeepVariant Pipeline

This pipeline performs an end-to-end variant calling, starting from raw fastQ files to a single-sample VCF. It has been pre-configured for the CCGA MedCluster. 

## Documentation

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Accuracy
Deepvariant has been shown to be highly accurate to the Genome-in-a-Bottle benchmark sets. Using in-house produced Illumina WGS data from GIAB sample NA12878, the following concordance scores were achieved within the GIAB high confidence intervals:


| Variant | Recall   | Precision |
| ------- | -------- | --------- |
| INDEL   | 0.990528 | 0.995649  |
| SNP     | 0.999513 | 0.999728  |
