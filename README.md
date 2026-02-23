# HCMV Transcriptome Pipeline Project

COMP 483 independent project comparing Human cytomegalovirus (HCMV) transcriptomes 2- and 6-days post-infection (dpi).

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [BioPython](https://biopython.org/wiki/Getting_Started)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [kallisto](https://pachterlab.github.io/kallisto/download)
- [R](https://www.r-project.org) and [sleuth](https://pachterlab.github.io/sleuth/download)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SPAdes](https://ablab.github.io/spades/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)
- [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)

All dependencies can also be installed by creating an environment in [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) with the needed tools/packages. See [Managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) and [sleuth installation via conda](https://github.com/pachterlab/sleuth)

## Running the Pipeline

### All Input Reads

For the complete pipeline (without test data), the rule `fasterq_dump` in the Snakefile runs the fasterq-dump command for each sample, creating fastq files that were used in the rest of the pipeline, resulting in `La_Rosa_PipelineReport.txt`.

This rule is still included in the Snakefile, but does not run since the rules that depend on fastq files as input will already have those files provided in `test_data/`

### Sample Input Reads

fastq files with only the first 10,000 reads per sample are provided in `test_data/`. After cloning the repo, run the full pipeline with this data using:

`snakemake -c1`

or

`snakemake -c2`

### Output

The pipeline will create a file called `PipelineReport.txt` in the root directory. This report contains the following information:
- Number of coding sequences (CDS) in the HCMV genome
- A tab-delimited table of the differentially expressed genes between the 2dpi and 6dpi samples
- The number of read pairs for each sample before and after Bowtie2 filtering
- A tab-delimited table showing the top 5 BLAST hits for each sample's longest contig from its assembly

### Cleanup
To remove all output files/folders created from running the pipeline:

`snakemake cleanup -c1`

Note that this will also remove the `PipelineReport.txt` generated from running the pipeline.

## References and Resources Used

- [Biopython SeqFeature documentation](https://biopython.org/docs/1.75/api/Bio.SeqFeature.html)
- [Biopython SeqFeature workshop](https://github.com/peterjc/biopython_workshop/blob/master/using_seqfeatures/README.rst)
- [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [Command line arguments in R scripts](https://wresch.github.io/2013/06/20/commandline-args-in-R.html#:~:text=args%20)
- [NCBI datasets HCMV genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000845245.1/)
- [Cheng et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/29158406)
- Lecture slides, example code, and other instructional materials courtesy of Dr. Heather Wheeler
