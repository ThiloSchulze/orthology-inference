#!/bin/bash

#readonly ARG_COUNT="$#" #see process_argument comment


usage() { printf "%s" "\
_____________________________________________________________________________________
There are four requirements in order to run this script:

  1) This script needs to be run in the same directory as your orthogroup files.
  2) Orthogroup sequence files need to be in fasta format and have the .fa extension.
  3) Fasta headers need to start with SPECIES_NAME@
  4) If you re-run the analysis in the same directory you need to provide -r
_____________________________________________________________________________________

usage:
  filtering_og.sh [--species] [--max_sequence]
example:
  ./filtering_og.sh -s 0.7 -m 300
description:
  Filter a set of orthogroup files for a) minimum % value of included species
                                       b) maximum number of included sequences
options:
  -h, --help          display this help message and exit
  -s, --species       determine a minimum % value of species in regard to max species diversity. (default: '$species')
                      only orthogroups with equal or higher % will be filtered. input value needs to be min 0 and max 1.
                      (e.g. with -s 0.75 and a dataset of 12 species, all orthogroups with 9 or more different species are filtered.)
  -m, --max_sequence  determine a maximum amount of sequences (default: '$sequence'). only orthogroups with equal to or less sequences will be filtered.
  -r, --rerun         delete existing fasta files in ./filtered_orthogroups/
"
        exit 1
        }

error() {
  local message="$1"
  local status="${2-1}" # default exit status: 1
  echo "filtering_og.sh error: $message
Please provide the following message to view the help menu:
./filtering_og.sh -h"
  exit "$status"
}


# Parse user-provided arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -h | --help)
      usage
      ;;
    -s | --species)
      species="$2"
      shift # past argument (both shifts are necessary to go from species to sequence)
      shift # past value
      ;;
    -m | --max_sequence)
      sequence="$2"
      shift # past argument
      shift # past value
      ;;
    -r | --rerun)
      rerun_status=1
      rm ./filtered_orthogroups/*.fa
      shift # past argument
      shift # past value
      ;;
    --) # end of all options
      break
      ;;
    -*) # unknown option
      error "unknown option: $key"
      ;;
    *) # end of options
      [[ -n "$basedir" ]] && error "unrecognized positional argument: $key"
      basedir="$key"
      shift # past argument
      ;;
  esac
done

#percentage_species=$( echo "$species ./* 100" | bc)
#echo "$percentage_species"
#process_argument() {
#  # Display an error message if no arguments were provided
#  [[ ARG_COUNT -lt 1 ]] && usage
#}
#process_argument    #Only un-comment this if you want that this program does not run without input

percentage_species=$( bc <<< "$species * 100" )

echo "
Sorting each orthogroup file into ./filtered_orthogroups/ that meets the following requirements:

1) Orthogroup includes at least ${percentage_species}% of species in the dataset.
2) Orthogroup includes at maximum $sequence total sequences.

Please view orthogroup_info_summary.txt or orthogroup_info_extended.txt for the results."

####Pre-filtering####

##Species diversity##
species_per_og=$( for file in *.fa
do 
  grep '^>' "$file" |\
  sed 's/@.*//' |\
  sort -n |\
  uniq |\
  wc -l
done )  #Searches each orthogroup file for its species diversity; saves it to the variable.

orthogroup_count=$( echo "$species_per_og" | wc -l ) #determines total orthogroup count
species_per_orthogroup=$( echo "$species_per_og" | sort -n | uniq -c ) #lists orthogroups separated by species diversity within

total_sum_species=$( echo "$species_per_og" | paste -sd+ - | bc )  #Calculating mean species diversity of orthogroups
mean_species_per_og=$( echo "scale=2 ; $total_sum_species / $orthogroup_count" | bc )
median_species_per_og=$( echo "$species_per_og" | sort -n | awk ' { a[i++]=$1; }
    END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' )

max_sp_div=$( echo "$species_per_og" | sort -nrk1,1 | head -1 | cut -d ' ' -f3 ) #determines highest species diversity out of all orthologs
cutoff=$( bc <<< "scale=3 ; $max_sp_div * $species" )  # percentage% of highest species diversity -> later cut off all orthogroups with lower species diversity
rounded_cutoff=$( echo "$cutoff" | awk '{print int($1+0.5)}' ) #Round species divesity cutoff accordingly

mkdir -p filtered_orthogroups   #-p means: only create folder if it doesn't already exist

##Sequences##

sequences_per_og=$( for file in *.fa
do 
  grep '^>' "$file" |\
  sort -n |\
  wc -l
done )  #Searches each orthogroup file for its species diversity; saves it to the variable.

sequences_per_orthogroup=$( echo "$sequences_per_og" | sort -n | uniq -c ) #lists orthogroups separated by total sequence count

total_sum_sequences=$( echo "$sequences_per_og" | paste -sd+ - | bc )  #Calculating mean sequence count of orthogroups
mean_sequences_per_og=$( echo "scale=2 ; $total_sum_sequences / $orthogroup_count" | bc )
median_sequences_per_og=$( echo "$sequences_per_og" | sort -n | awk ' { a[i++]=$1; }
    END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' )


for file in *.fa
do
  species_count=$( grep '^>' "$file" |\
  sed 's/@.*//' |\
  sort | uniq |\
  wc -l ) &&
  sequence_count=$( grep '^>' "$file" |\
  sort -n |\
  wc -l )
if (( species_count >= rounded_cutoff )) && (( sequence_count <= sequence ))
then cp "$file" "filtered_orthogroups/${file}"
fi
done #Copy all orthogroup files with species diversity above $species in % AND less than $sequence total sequences to ./filtered_orthogroups.

####Post-filtering####

##Species diversity##
species_per_og_pf=$( for file in ./filtered_orthogroups/*.fa
do 
  grep '^>' "$file" |\
  sed 's/@.*//' |\
  sort -n |\
  uniq |\
  wc -l
done )  #Searches each orthogroup file for its species diversity; saves it to the variable.

orthogroup_count_pf=$( echo "$species_per_og_pf" | wc -l ) #determines total orthogroup count
species_per_orthogroup_pf=$( echo "$species_per_og_pf" | sort -n | uniq -c ) #lists orthogroups separated by species diversity within

total_sum_species_pf=$( echo "$species_per_og_pf" | paste -sd+ - | bc )  #Calculating mean species diversity of orthogroups
mean_species_per_og_pf=$( echo "scale=2 ; $total_sum_species_pf / $orthogroup_count_pf" | bc )
median_species_per_og_pf=$( echo "$species_per_og_pf" | sort -n | awk ' { a[i++]=$1; }
    END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' )

##Sequences##

sequences_per_og_pf=$( for file in ./filtered_orthogroups/*.fa
do 
  grep '^>' "$file" |\
  sort -n |\
  wc -l
done )  #Searches each orthogroup file for its species diversity; saves it to the variable.

sequences_per_orthogroup_pf=$( echo "$sequences_per_og_pf" | sort -n | uniq -c ) #lists orthogroups separated by total sequence count

total_sum_sequences_pf=$( echo "$sequences_per_og_pf" | paste -sd+ - | bc )  #Calculating mean sequence count of orthogroups
mean_sequences_per_og_pf=$( echo "scale=2 ; $total_sum_sequences_pf / $orthogroup_count_pf" | bc )
median_sequences_per_og_pf=$( echo "$sequences_per_og_pf" | sort -n | awk ' { a[i++]=$1; }
    END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' )



echo "
Sorting each orthogroup file into ./filtered_orthogroups/ that meets the following requirements:

1) Orthogroup includes at least ${species}/1 of species in the dataset.
2) Orthogroup includes at maximum $sequence total sequences.

==================================
Stats on orthogroups pre-filtering
==================================

The total number of orthogroups is $orthogroup_count.


Stats on species diversity in orthogroups
#########################################

Mean species diversity per orthogroup is $mean_species_per_og.
Median species diversity per orthogroup is $median_species_per_og.

The distribution of orthogroups with X species is the following:
Number of Orthogroups | Number of different species within
$species_per_orthogroup


Stats on sequence abundance in orthogroups
##########################################

Mean sequence number per orthogroup is $mean_sequences_per_og.
Median sequence number per orthogroup is $median_sequences_per_og.

The distribution of orthogroups with X sequences is the following:
Number of Orthogroups | Number of sequences within

$sequences_per_orthogroup

 ___________________________________________________________________________________
|Sorting all orthogroups with: $rounded_cutoff or more different species            
|                              less than $sequence sequences into ./filtered_orthogroups/ 
|___________________________________________________________________________________

===================================
Stats on orthogroups post-filtering
===================================

The total number of orthogroups is $orthogroup_count_pf.


Stats on species diversity in orthogroups
#########################################

Mean species diversity per orthogroup is $mean_species_per_og_pf.
Median species diversity per orthogroup is $median_species_per_og_pf.

The distribution of orthogroups with X species is the following:
Number of Orthogroups | Number of different species within
$species_per_orthogroup_pf


Stats on sequence abundance in orthogroups
##########################################

Mean sequence number per orthogroup is $mean_sequences_per_og_pf.
Median sequence number per orthogroup is $median_sequences_per_og_pf.

The distribution of orthogroups with X sequences is the following:
Number of Orthogroups | Number of sequences within

$sequences_per_orthogroup_pf" > ./orthogroup_info_extended.txt


filtering_ogs_command="filtering_ogs --species ${species} --max_sequence $sequence"
  if [[ "$rerun_status" == 1 ]]
  then filtering_ogs_command+=" --rerun"
  fi

percentage_orthogroup=$( bc -l <<<  "scale=2; ($orthogroup_count_pf / $orthogroup_count ) * 100" )

echo "
Input
=====
$filtering_ogs_command

Output
======
Filtered Orthogroups include at least ${rounded_cutoff} different species and have $sequence sequences at maximum.
The total number of orthogroups is $orthogroup_count (pre-filtering) and $orthogroup_count_pf (post-filtering). ${percentage_orthogroup}% of the dataset was filtered.

Mean species diversity per orthogroup is $mean_species_per_og (pre) and $mean_species_per_og_pf (post).
Median species diversity per orthogroup is $median_species_per_og (pre) and $median_species_per_og_pf (post).

Mean sequence number per orthogroup is $mean_sequences_per_og (pre) and $mean_sequences_per_og_pf (post).
Median sequence number per orthogroup is $median_sequences_per_og (pre) and $median_sequences_per_og_pf (post).

Filtered orthologs have been sorted into ./filtered_orthogroups/" > ./orthogroup_info_summary.txt