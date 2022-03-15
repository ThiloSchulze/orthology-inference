#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

// Direct user to help message
if ( ! params.dataset ) exit 0, guideMessage()

// Creation of dataset directory channel with file extension filter and check if empty
ch_sequences = Channel
  .fromPath( params.dataset, type: 'dir' )
  //.filter { it =~/.*\.fa|.*\.faa|.*\.fas|.*\.pep/ }  //Does not work
  .view()
  .ifEmpty { exit 1,
             "No FASTA files found with pattern '${params.dataset}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fa', '.faa', '.fasta', '.fas', or '.pep'\n" +
             "Please provide a directory that includes your dataset." }

def helpMessage() {
    log.info"""
    Basic usage example:
        (1) Amino acid data:
        nextflow run main.nf --dataset 'directory/' --m LG
        (2) Nucleotide data:
        nextflow run main.nf --dataset 'directory/' --d --m HKY+F

    Application example:
        nextflow run main.nf --dataset 'directory/' --filter_species 0.8 --max_sequence 300 --maxiterate 1000 --ufboot 1000 --m MFP

    Description:
        A pipeline for graph-based and tree-based orthology inference of assembled WGS data

    Pipeline summary:
        1. Graph-based orthology inference using OrthoFinder
        2. Removal of orthogroups with low information content
        3. Alignment of gene trees using mafft
        4. Gene tree inference using IQTREE
  
    Mandatory arguments:
        --dataset           path to directory containing your input fasta files. One file per species.
                            (valid file types: '.fa', '.faa', '.fasta', '.fas', or '.pep')

    Input/output options:
        --output            specifies path the output is saved to (default: $params.output)

    Resource allocation:
        --memory            memory used by the pipeline  (default: $params.memory)
        --threads           number of threads to utilize (default: $params.threads)
        --time          runtime of the pipeline      (default: $params.time)

    Graph-based orthology inference (OrthoFinder):
        --d                 Input is DNA sequences
        --I NUM             MCL inflation parameter (default: 1.5)
        --S                 Sequence search program (default: diamond)
                            (options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl)

    Filtering orthogroups (filtering_ogs):
        --filter_species    Remove orthogroups including less than X% of of species in dataset (default: $params.filter_species)
        --max_sequence NUM  Maximum number of sequences per orthogroup (default: $params.max_sequence)

    Aligning gene trees (Mafft):
        --6merpair          Distance is calculated based on the number of shared 6mers (default: on)
        --globalpair        Pairwise alignments are computed with the Needleman-Wunsch algorithm
        --localpair         Pairwise alignments are computed with the Smith-Waterman algorithm 
        --genafpair         Pairwise alignments are computed with a local algorithm with the generalized affine gap cost (Altschul 1998)
        --fastapair         Pairwise alignments are computed with FASTA (Pearson and Lipman 1988)
        --memsave           Use the Myers-Miller (1988) algorithm
        --weighti NUM       Weighting factor calculated from pairwise alignments (default: 2.7)
                            (valid with globalpair, localpair, genafpair, fastapair)
        --retree NUM        Guide tree is built number times in the progressive stage (valid with 6mer distance, default: 2)
        --maxiterate NUM    X cycles of iterative refinement are performed (default: 0)
        --nofft             Do not use FFT approximation in group-to-group alignment. (default: off)
        --parttree          Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the 6mer distance (default: off)
                            (options: parttree, dpparttree, fastaparttree)
        --quiet_mafft       Do not report progress
        --op NUM            Gap opening penalty (default: 1.53)
        --op NUM            Offset (works like gap extension penalty, default: 0.0)

    Infering gene trees (IQTREE):
        --m                 Choose specific model or choose model selection [e.g. TEST, MF, MFP] (default: $params.m)
        --quiet_iqtree      Do not report progress
        --epsilon NUM       Likelihood epsilon for parameter estimate (default 0.01)
        --ninit NUM         Number of initial parsimony trees (default: 100)
        --ntop NUM          Number of top initial trees (default: 20)
        --nbest NUM         Number of best trees retained during search (defaut: 5)
        --ufboot NUM        Replicates for ultrafast bootstrap (>=1000)
        --ufjack NUM        Replicates for ultrafast jackknife (>=1000)
        --jackprop NUM      Subsampling proportion for jackknife (default: 0.5)
        --nmax NUM          Maximum number of iterations (default: 1000)
        --nstep NUM         Iterations for UFBoot stopping rule (default: 100)
        --bcor NUM          Minimum correlation coefficient (default: 0.99)
        --beps NUM          RELL epsilon to break tie (default: 0.5)
        --bnni              Optimize UFBoot trees by NNI on bootstrap alignment
        --boot NUM          Replicates for bootstrap + ML tree + consensus tree
        --jack NUM          Replicates for jackknife + ML tree + consensus tree
        --alrt NUM          Replicates for SH approximate likelihood ratio test (e.g. 1000)
        --alrt0             Parametric aLRT test (Anisimova and Gascuel 2006)
        --abayes            approximate Bayes test (Anisimova et al. 2011)
        --lbp NUM           Replicates for fast local bootstrap probabilities

    Miscellaneous:
        --help              display this help message and exit
        --version           display the pipeline's version number and exit
    """.stripIndent()
}

def guideMessage() {
    log.info"""
    No dataset was provided. Please type the following message for more information:
    nextflow run main.nf --help
    """.stripIndent()
}

process orthofinder {
    // Copies output file in output folder
    publishDir "${params.output}/OrthoFinder", mode: 'copy'

    input:
    // All files in the specified --dataset 'directory/'
    path(species)

    output:
    // Direct each orthogroup file to filtering
    path "${params.dataset}OrthoFinder/Results_$params.orthofinder_output_dir/*" optional true
    path "${params.dataset}OrthoFinder/Results_$params.orthofinder_output_dir/Orthogroup_Sequences/*", emit: 'ogs'
    path "${params.dataset}OrthoFinder/Results_$params.orthofinder_output_dir/OrthoFinder_command.txt", emit: command

    script:
    // Align sequences within each orthogroup file

    flagsOrthoFinder = "-f $species -M msa -A mafft -os -n $params.orthofinder_output_dir -t $task.cpus -a $params.a"
    if ( params.d )
        flagsOrthoFinder += " -d"
    if ( params.I )
        flagsOrthoFinder += " -I $params.I"
    if ( params.S )
        flagsOrthoFinder += " -S $params.S"
    commandOrthoFinder = "orthofinder.py $flagsOrthoFinder"

    """
    $commandOrthoFinder
    echo "$commandOrthoFinder" > '${params.dataset}OrthoFinder/Results_$params.orthofinder_output_dir/OrthoFinder_command.txt'
    """

}

process filtering {
    publishDir "${params.output}/filtering", mode: 'copy'

    beforeScript 'ulimit -s unlimited'

    input:
    // All orthogroup sequences derived from the selected dataset.
    path(orthogroups)

    output:
    // Direct filtered orthogroup files to mafft
    path "filtered_orthogroups/*", emit: filtered_ogs
    path "*.txt"
    path "*_command.txt", emit: command

    script:
    // Remove orthogroups that are uninformative due to a) a low number of different species or b) too many sequences



    """
    filter_ogs.sh --species $params.filter_species --max_sequence $params.max_sequence --rerun
    echo "filter_ogs.sh --species $params.filter_species --max_sequence $params.max_sequence --rerun" > filter_ogs_command.txt
    """

}


process mafft {
    publishDir "${params.output}/mafft", mode: 'copy'

    input:
    // All filtered orthogroups
    path(filtered_orthogroups)

    output:
    // Direct orthogroups with aligned sequences to iqtree and phylopypruner
    path "*_mafft.fa", emit: ogs_aligned
    path "*.txt"
    path "*_command.txt", emit: command
    path "list_of_species.txt", emit: species_list

    script:
    // Use mafft to align sequences within each orthogroup file
  flagsmafft = "--thread $task.cpus"
  if ( params.globalpair )
    flagsmafft += " --globalpair"
  if ( params.localpair )
    flagsmafft += " --localpair"
  if ( params.genafpair )
    flagsmafft += " --genafpair"
  if ( params.fastapair )
    flagsmafft += " --fastapair"
  if ( params.memsave )
    flagsmafft += " --memsave"
  if ( params.weighti  )
    flagsmafft += " --weighti $params.weighti"
  if ( params.retree )
    flagsmafft += " --retree $params.retree"
  if ( params.maxiterate )
    flagsmafft += " --maxiterate $params.maxiterate"
  if ( params.nofft )
    flagsmafft += " --nofft"
  if ( params.parttree )
    flagsmafft += " --parttree"
  if ( params.dpparttree )
    flagsmafft += " --dpparttree"
  if ( params.fastaparttree )
    flagsmafft += " --fastaparttree"
  if ( params.quiet_mafft )
    flagsmafft += " --quiet"
  if ( params.op )
    flagsmafft += " --op $params.op"
  if ( params.ep )
    flagsmafft += " --ep $params.ep"
  commandmafft = "$params.mafft $flagsmafft"

    """
    for fa in ${filtered_orthogroups}; do $commandmafft "\$fa" > \${fa/.*/_mafft.fa}; done 
    echo "$commandmafft" > 'mafft_command.txt'
    grep --no-filename '^>' *.fa | tr -d '>' | sed 's/@.*//g' | sort | uniq > list_of_species.txt
    """
}

//  for faa in ${og}; do
//  mafft \$faa > \${faa/.faa/_mafft.faa}
//  done
//"phylopypruner_prep/\${fa/_mafft.fa/.fa}"
//    """
//    for fa in ${filtered_orthogroups}; do mafft --thread $task.cpus "\$fa" > \${fa/.*/_mafft.fa}; done 
//    echo "for fa in ${filtered_orthogroups}; do mafft "\$fa" > \${fa/.*/_mafft.fa}; done" > mafft_command.txt
//    """

process iqtree {
    publishDir "${params.output}/iqtree", mode: 'copy'

    input:
    // All filtered orthogroup files with aligned sequences
    path(og_alignment)

    output:
    // Direct infered gene trees to phylopypruner
    path "*"
    path "*.treefile", emit: gene_tree_files
    path "*_command.txt", emit: command

    script:
    // Infer a gene tree from each filtered orthogroup file
  flagsiqtree = "-T $task.cpus"
  if ( params.quiet_iqtree )
    flagsiqtree += " --quiet"
  if ( params.epsilon )
    flagsiqtree += " --epsilon $params.epsilon"
  if ( params.ninit )
    flagsiqtree += " --ninit $params.ninit"
  if ( params.ntop )
    flagsiqtree += " --ntop $params.ntop"
  if ( params.nbest )
    flagsiqtree += " --nbest $params.nbest"
  if ( params.ufboot )
    flagsiqtree += " --ufboot $params.ufboot"  // >=1000
  if ( params.ufjack )
    flagsiqtree += " --ufjack $params.ufjack"  // >=1000
  if ( params.jackprop )
    flagsiqtree += " --jack-prop $params.jackprop" // 0.5
  if ( params.nmax )
    flagsiqtree += " --nmax $params.nmax"  // 1000
  if ( params.nstep )
    flagsiqtree += " --nstep $params.nstep" // 100
  if ( params.bcor )
    flagsiqtree += " --bcor $params.bcor"  // 0.99
  if ( params.beps )
    flagsiqtree += " --beps $params.beps"  // 0.5
  if ( params.bnni )
    flagsiqtree += " --bnni"
  if ( params.boot )
    flagsiqtree += " --boot $params.boot"  // >=1000
  if ( params.jack )
    flagsiqtree += " --jack $params.jack"  // >=1000
  if ( params.alrt )
    flagsiqtree += " --alrt $params.alrt"  // >=1000
  if ( params.alrt0 )
    flagsiqtree += " --alrt 0"
  if ( params.abayes )
    flagsiqtree += " --abayes"
  if ( params.lbp )
    flagsiqtree += " --lbp $params.lbp"
  if ( params.m )
    flagsiqtree += " -m $params.m"
  commandiqtree = "iqtree $flagsiqtree"

    """
    for alignment in ${og_alignment}; do $commandiqtree -s "\$alignment"; done 
    echo "$commandiqtree" > iqtree_command.txt
    """

}

process iq_to_ppp {
    publishDir "${params.output}/iq_to_ppp", mode: 'copy'

input:
    // All orthogroup sequences derived from the selected dataset.
    path(species_list)
    path(tre)
    path(mafft_alignments)

    output:
    // Direct filtered orthogroup files to mafft
    path "phylopypruner_prep", type: 'dir'
    path "*"


    script:
    // Remove orthogroups that are uninformative due to a) a low number of different species or b) too many sequences

    """
    mkdir -p phylopypruner_prep
    for filename in $tre
    do
      while read -r otu
      do
        sed -i "s/\${otu}_/\${otu}@/g" "\${filename}"
      done < "$species_list"
      mv "\$filename" "phylopypruner_prep/\${filename/_mafft.fa.treefile/.tre}"
    done

    for mafft_fa in $mafft_alignments
    do mv "\$mafft_fa" "phylopypruner_prep/\${mafft_fa/_mafft.fa/.fa}"
    done

    phylopypruner --dir 'phylopypruner_prep/'
    """
//     phylopypruner --dir 'phylopypruner_prep/'
}


process fasttrees {
    publishDir "${params.output}/FastTree", mode: 'copy'

    input:
    // All filtered orthogroup files with aligned sequences
    path(og_alignment)

    output:
    // Direct infered gene trees to phylopypruner
    path "*.tre", emit: tre_files
    path "*_command.txt", emit: command

    script:


    """
    for alignment in ${og_alignment}; do fasttree < "\$alignment" > "\${alignment%_mafft*}.tre"; done 
    echo "fasttree < "\$alignment" > \${alignment%_mafft*}.tre" > FastTree_command.txt
    """

}


process fast_to_ppp {
    publishDir "${params.output}/fast_to_ppp", mode: 'copy'

    input:
    // All orthogroup sequences derived from the selected dataset.
    path(tre_file)
    path(alignments)

    output:
    // Direct filtered orthogroup files to mafft
    path "phylopypruner_prep", type: 'dir', emit: phylo_prep
    path "*"


    script:
    // Remove orthogroups that are uninformative due to a) a low number of different species or b) too many sequences

    """
    mkdir -p phylopypruner_prep
    for tre in $tre_file
    do
      mv "\$tre" "phylopypruner_prep/\$tre"
    done

    for alignment in $alignments
    do 
      mv "\$alignment" "phylopypruner_prep/\${alignment/%_mafft.fa}.fa"
    done
    """
  //  phylopypruner --dir 'phylopypruner_prep/'
}

process ppp {
    publishDir "${params.output}/phylopypruner", mode: 'copy'

    input:
    // All orthogroup sequences derived from the selected dataset.
    path(pruner_prep)

    output:
    path "phylopypruner_output", type: 'dir'
    path "*_command.txt", emit: command


    script:
    // Infer a gene tree from each filtered orthogroup file
  flagsphylopypruner = "--dir $pruner_prep --threads $params.threads --output phylopypruner_output"
  if ( params.no_plot )
    flagsphylopypruner += " --no-plot"
  if ( params.no_supermatrix )
    flagsphylopypruner += " --no-supermatrix"
  if ( params.min_len )
    flagsphylopypruner += " --min-len $params.min_len"
  if ( params.trim_lb )
    flagsphylopypruner += " --trim-lb $params.trim_lb"
  if ( params.min_pdist )
    flagsphylopypruner += " --min-pdist $params.min_pdist"
  if ( params.min_support )
    flagsphylopypruner += " --min-support $params.min_support"
  if ( params.trim_divergent )
    flagsphylopypruner += " --trim-divergent $params.trim_divergent"
  if ( params.trim_freq_paralogs )
    flagsphylopypruner += " --trim-freq-paralogs $params.trim_freq_paralogs"
  if ( params.include )
    flagsphylopypruner += " --include $params.include"
  if ( params.outgroup )
    flagsphylopypruner += " --outgroup $params.outgroup"
  if ( params.root )  
    flagsphylopypruner += " --root $params.root"
  if ( params.prune )
    flagsphylopypruner += " --prune $params.prune" //LS
  if ( params.force_inclusion )
    flagsphylopypruner += " --force-inclusion $params.force_inclusion"
  if ( params.min_taxa )
    flagsphylopypruner += " --min-taxa $params.min_taxa"
  if ( params.min_otu_occupancy )
    flagsphylopypruner += " --min-otu-occupancy $params.min_otu_occupancy"
  if ( params.min_gene_occupancy )
    flagsphylopypruner += " --min-gene-occupancy $params.min_gene_occupancy"
  if ( params.subclades )
    flagsphylopypruner += " --subclades $params.min_otu_occupancy"
  if ( params.jackknife )
    flagsphylopypruner += " --jackknife $params.jackknife"
  commandphylopypruner = "phylopypruner $flagsphylopypruner"

    """
    $commandphylopypruner
    echo "$commandphylopypruner" > phylopypruner_command.txt
    """
  //  phylopypruner --dir 'phylopypruner_prep/'
}



process settings {
    publishDir "${params.output}/commands", mode: 'copy'

    input:
    path(command1)
    path(command2)
    path(command3)
    path(command4)
    path(command5)

    output:
    // Complete list of pipeline settings
    path "*.txt"

    script:
    """
    command_list=\$( cat $command1 $command2 $command3 $command4 $command5 )
    echo "Input
    =====

    \$command_list" > command_list.txt
    """

}

workflow {
    orthofinder(ch_sequences)
    filtering(orthofinder.out.ogs.collect())
    mafft(filtering.out.filtered_ogs)
    if (params.fasttree) {
      fasttrees(mafft.out.ogs_aligned)
      fast_to_ppp(fasttrees.out.tre_files, mafft.out.ogs_aligned)
      ppp(fast_to_ppp.out.phylo_prep)
      settings(orthofinder.out.command, filtering.out.command, mafft.out.command, fasttrees.out.command, ppp.out.command)
    }
    else {
      iqtree(mafft.out.ogs_aligned)
      iq_to_ppp(mafft.out.species_list, iqtree.out.gene_tree_files, mafft.out.ogs_aligned)
      ppp(fast_to_ppp.out.phylo_prep)
      settings(orthofinder.out.command, filtering.out.command, mafft.out.command, iqtree.out.command)
    }

//    settings(orthofinder.out.command, filtering.out.command, mafft.out.command, iqtree.out.command)
}
