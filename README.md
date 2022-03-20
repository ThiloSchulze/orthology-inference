> orthology-inference is still in development.

# Orthology Inference Pipeline
This pipeline performs graph-based and tree-based orthology inference of a dataset containing sequences of at least four species. Steps:

1. Graph-based orthology inference ([OrthoFinder](https://github.com/davidemms/OrthoFinder))
2. Filtering of uninformative orthogroups (custom script)
3. Multiple sequence alignment ([MAFFT](https://mafft.cbrc.jp/alignment/software/))
4. Inference of gene trees ([FastTree](http://www.microbesonline.org/fasttree/))
5. Tree-based orthology inference ([PhyloPyPruner](https://pypi.org/project/phylopypruner/))


## Installing nextflow

[Nextflow](https://www.nextflow.io/) is required to run this pipeline. For
installing Nextflow, simply follow the [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)
in their documentation.

## Downloading this pipeline

To download this pipeline, go to the folder where you wish to store this
pipeline and then run the following:

    $ git clone https://github.com/ThiloSchulze/orthology-inference

## Preparations before usage

1. Installation of Nextflow and the programs listed as steps.
2. Input for the pipeline must be a single directory containing your input fasta files, with one file per species.

   The following file extensions are recognized: `.fa`, `.faa`, `.fasta`, `.fas` and  `.pep`.
3. Each header within a fasta file must start by an identifier unique to the species and be followed by @, e.g. `>GENUS_SPECIES@`.

## Basic parameters

#### `--dataset`

* Mandatory argument
* Specifies the location of the directory containing your fasta files.
* Example: `--dataset 'directory/'`

#### `--output`

* Specifies output directory for generated results.
* If the specified folder is not existant, it will be created.
* By default it is set to: `inference_out`

#### `--filter_species`

* Choose a percentage value as decimal. Orthogroups containing less species than specified will be exempt from the analysis.
* Example: Your dataset contains 20 species. `--filter_species 0.75` will remove all orthogroups containing markers from less than 15 species.
* Aims to reduce runtime of MAFFT and FastTree
* By default it is set to: `--filter_species 1`

#### `--max_sequence`

* Choose a number. Orthogroups containing equal or more sequences will be exempt from the analysis.
* Example: `--max_sequence 150` will remove all orthogroups containing 150 or more sequences.
* Aims to remove Orthogroups containing sequencing errors.
* By default it is set to: `--max_sequence 200`

#### `--help`

* Executes the help menu
* Example: `nextflow run main.nf --help`


#### Settings for MAFFT, FastTree / IQTREE, PhyloPyPruner

This pipeline is meant to be customizable. Most options available for each of the listed programs are also available in this pipeline.
The input options were kept as similar to the respective program as possible. If you wish to create the gene trees with IQTREE instead of FastTree,
please provide `--iqtree` (careful: the runtime will be increased drastically).

For specific parameter options please set the `--help` parameter in the base directory of the pipeline.

## Usage

#### Basic usage example:

(1) Amino acid data:
`nextflow run main.nf --dataset 'directory/'`
    
(2) Nucleotide data:
`nextflow run main.nf --dataset 'directory/' --d`


#### Application example:

`nextflow run main.nf --dataset 'directory/' --filter_species 0.8 --max_sequence 300 --maxiterate 1000 --min_len 100 --trim_lb 5 --min_support 0.70 --min_taxa 10 --min_otu_occupancy 0.1 --min_gene_occupancy 0.1 --threads 4`
