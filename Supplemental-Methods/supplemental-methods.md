# Supplemental Code

**Authors: Jonah Tang and Jonathon L. Baker**

This document catalogues the scripts and computational methods used in "The salivary virome during childhood dental caries". Saliva samples were collected from individuals as detailed in a previous study by our group, Baker et al. 2021 (PMID: 33239396). While the previous study examined metagenomics of these samples with a primarily bacterial focus, this study examined the oral virome and its relationship to caries, the bacteriome, and host immunological markers. This study also identified novel viral taxa. The raw sequencing reads of the oral metagenomes and are available on NCBI with accession numbers PRJNA478018 and SRP151559. The bacterial genomes assembled from those reads  are available on NCBI under the accession number PRJNA624185.  The 35 vMAGs representing novel vOTUs rated as either ‘complete’ or ‘high-quality’ by CheckV are available on NCBI with the accession number PRJNA478018.

## Abstract
While many studies have examined the bacterial taxa associated with dental caries, the most common chronic infectious disease globally, little is known about the caries-associated virome. In this study, the salivary viromes of 21 children with severe caries (>2 dentin lesions) and 23 children with healthy dentition were examined. 2,485 viral metagenome-assembled genomes (vMAGs) were identified, binned, and quantified from the metagenomic assemblies. These vMAGs were mostly phage, and represented 1,865 unique species-level viral operational taxonomic units (vOTUs), 478 of which appear to be novel.  The metagenomes were also queried for all 3,858 unique species-level vOTUs of DNA viruses with a human host on NCBI Virus, however all but *Human betaherpesvirus 7* were at very low abundance in the saliva. The oral viromes of the children with caries exhibited significantly different beta diversity compared to the oral virome of the children with healthy dentition; several vOTUs predicted to infect *Haemophilus* and *Neisseria* were strongly correlated with health, and a few vOTUs predicted to infect Saccharibacteria and *Veillonella* were correlated with caries. Co-occurrence analysis indicated that phage typically co-occurred with both their predicted hosts and with bacteria that were themselves associated with the same disease status. Overall, this study provided the sequences of 35 complete or nearly complete novel oral phages and illustrated the potential significance of the oral virome in the context of dental caries, which has been largely overlooked. This work represents an important step towards the identification and study of phage therapy candidates which treat or prevent caries pathogenesis.

## ViWrap

Metagenomic assemblies were first processed with [ViWrap](https://github.com/AnantharamanLab/ViWrap) v1.3.0, a comprehensive tool for identification, binning, quality analysis, taxonomic prediction, and host prediction of vMAGs. Each sample was processed with the following arguments:

```sh
conda run -p ViWrap_conda_environments/ViWrap ViWrap/ViWrap run \
            --input_metagenome virome/mspades${sample_id}-contigs.fasta \
            --input_reads "$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_1.fastq),$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_1.fastq)" \
            --out_dir virome_output/slurm-${SLURM_ARRAY_JOB_ID}/${sample_id}_metaG-ViWrap_out \
            --db_dir miniforge3/envs/viwrap_env/ViWrap/ViWrap_db \
            --identify_method vb-vs \
            --conda_env_dir ViWrap_conda_environments/ \
            --threads 24 \
            --input_length_limit 5000
```

The output provided binning data from the vRhyme component of the ViWrap pipeline. However, this tool did not yield output for sample 23, resulting in its absence from certain downstream analyses. While it did not yield any binned vMAGs, sample 23 was still used for the read mapping steps using BWA-MEM and CoverM described later in this article.

### Compilation of ViWrap Output

From the ViWrap output directory `08_ViWrap_summary_outdir` for each sample, the following files of tabulated data were joined together into one file (`<sampleid>_Compiled_Output.txt`):

- `Host_prediction_to_genome_m90.csv`
- `Host_prediction_to_genus_m90.csv`
- `Tax_classification_result.txt`
- `Virus_normalized_abundance.txt`
- `Virus_raw_abundance.txt`
- `Species_cluster_info.txt`

This compilation was performed with the following script:

#### viwrap_compile_tables.sh

```sh
#!/bin/bash

# Performs a full outer join on tabulated ViWrap outputs.

# Parse arguments
sample_id=$1
if [ "$1" -eq 22 ] || [ "$1" -eq 31 ]; then
    sample_id="${1}A"
fi

# File paths
main_path="virome/${sample_id}_metaG-ViWrap_out"
summary_data_file="${main_path}/08_ViWrap_summary_outdir/Virus_summary_info.txt"
output_file="ViWrap_Compiled_Output/${sample_id}_Compiled_Output.txt"
tmp_output=$(mktemp)
tmp_header=$(mktemp)
tmp_tsv=$(mktemp)

# Converts a csv file to tsv, saving it the tmp_tsv variable
# Arguments:
#   1: Path to file to convert
convert_to_tsv () {
    tr ',' '\t' < $1 > $tmp_tsv
}

# Takes a tsv file and adds a header with the appropriate number of columns, saving it to the tmp_header variable
# Arguments:
#   1: Path to file
add_header () {
    local file_name=$(basename "$1")
    local header_line=""
    for i in $( seq 1 $(grep -o $'\t' <(tail -n 1 "${1}") | wc -l));
    do
        header_line+=$'\t'
        header_line+=$(echo "${file_name%.*}")_${i}
    done
    echo -e "${header_line}\n$(cat $1)" > $tmp_header
}

# Takes two files and joins them, writing to stdout
# Files must be in tsv format with a header line. Data must be sorted
# Arguments:
#   1: Path to first file
#   2: Path to second file
full_outer_join () {
    join --header -a 1 -a 2 -t $'\t' -o 'auto' -j 1 <(head -n 1 $1 && tail -n +2 $1 | sort --stable -k1,1) <(head -n 1 $2 && tail -n +2 $2 | sort --stable -k1,1)
}

# Create a temporary output file as a copy of the summary file
cp $summary_data_file $tmp_output

# Join contents of "Host_prediction_to_genome_m90.csv"
convert_to_tsv "${main_path}/08_ViWrap_summary_outdir/Host_prediction_to_genome_m90.csv"
full_outer_join $tmp_output $tmp_tsv | tee $output_file
cp $output_file $tmp_output

# Join contents of "Host_prediction_to_genus_m90.csv"
convert_to_tsv "${main_path}/08_ViWrap_summary_outdir/Host_prediction_to_genus_m90.csv"
full_outer_join $tmp_output $tmp_tsv | tee $output_file
cp $output_file $tmp_output

# Join contents of "Tax_classification_result.txt"
add_header "${main_path}/08_ViWrap_summary_outdir/Tax_classification_result.txt"
full_outer_join $tmp_output $tmp_header | tee $output_file
cp $output_file $tmp_output

# Join contents of "Virus_normalized_abundance.txt"
full_outer_join $tmp_output "${main_path}/08_ViWrap_summary_outdir/Virus_normalized_abundance.txt" | tee $output_file
cp $output_file $tmp_output

# Join contents of "Virus_raw_abundance.txt"
full_outer_join $tmp_output "${main_path}/08_ViWrap_summary_outdir/Virus_raw_abundance.txt" | tee $output_file
cp $output_file $tmp_output

# Join contents of "Species_cluster_info.txt"
convert_to_tsv "${main_path}/08_ViWrap_summary_outdir/Species_cluster_info.txt"
full_outer_join $tmp_output $tmp_tsv | tee $output_file
cp $output_file $tmp_output

# Clean up temporary files
rm $tmp_output
rm $tmp_tsv
rm $tmp_header
```

`viwrap_compile_tables.sh` was ran using the following script:

#### viwrap_compile_tables_run.sh

```sh
#!/bin/bash

for i in $(seq 1 49)
do
    if [ $i -ne 2 ] && [ $i -ne 23 ] && [ $i -ne 24 ] && [ $i -ne 32 ] && [ $i -ne 45 ]; then
        bash "viwrap_compile_tables.sh" "$i"
    fi
done
```

Finally, each `<sample-id>_Compiled_Output.txt` was compiled into a singular file `All_Compiled_Output.txt` with the following script; importantly, this script prepends the sample identifier to the names of the bins, providing an intermediate name for the genome (e.g. sample 1's `vRhyme_bin_1` becomes `1_vRhyme_bin_1`):

#### viwrap_compile_all_tables.sh

```sh
#!/bin/bash

# Compiles sets of data files in the same directory with matching headers, prepending the data rows with the appropriate identifier

# Prepends a string to each line of stdin, redirecting to stdout
# Arguments:
#   1: String to prepend to line
prepend () {
    while read line; do
        echo "${1}${line}"
    done
}

# Create temporary files
tmp_all_output=$(mktemp)

# Iterate through the files, copying the header line from the first file
# This works because we expect all data to align to the same header
first_loop=true
for file in ViWrap_Compiled_Output/*.txt; do
    if [[ "$first_loop" == true ]]; then
        head -n 1 $file | tee $tmp_all_output
        first_loop=false
    fi
    # Prepend the prefix of the file's name to the first column of each row of data
    file_prefix=$(basename $file .txt | cut -d_ -f1)
    tail -n +2 $file | prepend "${file_prefix}_" >> $tmp_all_output
done

# Copy over temporary file to output file
cp $tmp_all_output "ViWrap_Compiled_Output/All_Compiled_Output.txt"

# Clean up temporary files
rm $tmp_all_output
```

## Dereplication & Clustering

The assembled genomes were then dereplicated and organized into clusters using a 95% cutoff for average nucleotide identity (ANI) and 85% alignment fraction (AF), as specified by the MIUViG standards (Roux et al. 2019; PMID:30556814). We then selected a genome from each cluster to represent the species as a viral operational taxonomic unit (vOTU) for downstream analyses; selection was based on genome completeness as predicted by CheckV.

### Dereplication

Genomes were dereplicated into clusters using the [anvi'o](https://anvio.org) (development version) [`anvi-dereplicate-genomes`](https://anvio.org/help/7/programs/anvi-dereplicate-genomes/) program with the following arguments:

```sh
conda run -p "/miniforge3/envs/anvio-dev" anvi-dereplicate-genomes \
    -f "/derep/derep_virus_input.txt" \
    -o "/derep/derep_output_fastani/anvio_genome_similarity" \
    --program fastANI \
    --similarity-threshold 0.95 \
    --min-alignment-fraction 0.85 \
    --log-file "/derep/derep_output_fastani/anvio_derep_log.txt" \
    --num-threads 24 \
    --force-overwrite
```

It should be noted that the above input parameter (the `-f` option) takes a [fasta-txt artifact](https://anvio.org/help/7/artifacts/fasta-txt/), which we generated from a directory of FASTA files using the following script:

```sh
#!/bin/bash

output_file=derep_fasta_input.txt

rm $output_file
touch $output_file
echo -e "name\tpath" >> $output_file

find viwrap_fasta/*.fasta | while read line; do
    sample_id=$(basename $line .fasta)
    echo -e "${sample_id}\t${line}" >> $output_file
done
```

The `anvi-dereplicate-genomes` program provided a `CLUSTER_REPORT.txt` output file, containing the genomes that comprise each cluster.

### Clustering

Although `anvi-dereplicate-genomes` selects a cluster-representative genome based on centrality, the ViWrap output data indicated that some of these genomes had low percentages of completeness. To avoid a compounding loss of quality when mapping reads to the cluster-representative genomes when creating our OTU table, we opted to select for new cluster-representative genomes based on genome completeness percentage.

This process was performed with a script, where cluster-representative genomes are selected on attributes of the following priority:

| Priority | Category | `All_Compiled_Output.txt` Column Name |
| --- | --- | --- |
| 1 | Completeness Percentage | `completeness` |
| 2 | Confidence Level | `completeness_method` |
| 3 | Genome Size | `genome_size` |

For example, if two genomes in a cluster have identical values for completeness percentage, their confidence levels are checked against each other; the one with a higher level of confidence is selected as the cluster's representative genome.

#### find_most_complete_genome_in_cluster.sh

```sh
#!/bin/bash

# Takes anvio's CLUSTER_REPORT.txt file and our ViWrap All_Compiled_Data.txt file
# to find the best representative genome of the cluster by these attributes: Completeness > Confidence > Genome Size
# Arguments:
#   1: Anvio cluster data ("CLUSTER_REPORT.txt")
#   2: Compiled ViWrap data ("All_Compiled_Output.txt")
#   3: Output file path (including file name)

echo "Initializing files..."

# Create temporary files
# Tab delimited files' tabs are converted to unit separators to allow for proper field separation
# (Otherwise, rows with empty columns are improperly iterated through)
tmp_dir=$(mktemp -d)
tmp_cluster_report_data_translated=$(mktemp -p $tmp_dir)
tr '\t' '\037' < $1 > $tmp_cluster_report_data_translated
tmp_viwrap_data_translated=$(mktemp -p $tmp_dir)
tr '\t' '\037' < $2 > $tmp_viwrap_data_translated
# Create output file with and append header line
tmp_output=$(mktemp -p $tmp_dir)
echo -e "cluster_id\tcluster_size\tcluster_rep\tcluster_rep_completeness_value\tcluster_rep_confidence\tcluster_rep_confidence_method\tcluster_rep_genome_size\tcluster_rep_condition" >> $tmp_output
# Track temporary files for deletion on script exit
trap 'rm -rf -- "$tmp_dir"' EXIT

# Initialize array for determining confidence
declare -A enum_confidence_array=( ["(lower-bound)"]=1 ["(medium-confidence)"]=2 ["(high-confidence)"]=3 )

# Set Internal Field Separator to unit separator character
whitespace_ifs=$IFS
IFS=$'\037'

# Number of clusters (for debugging)
let total_clusters=$(sed -n '$=' $tmp_cluster_report_data_translated)-1

echo "Initializing arrays..."

# Convert files to arrays
readarray cluster_report_data < <(tail -n +2 $tmp_cluster_report_data_translated)
readarray viwrap_data < <(tail -n +2 $tmp_viwrap_data_translated)

# Create associative arrays for fast lookup
declare -A genome_size_aa
declare -A genome_completeness_value_aa
declare -A genome_completeness_confidence_method_aa
declare -A genome_completeness_confidence_value_aa

# Populate associative arrays from ViWrap data
for i in "${!viwrap_data[@]}"; do
    read -a viwrap_row_array < <(echo "${viwrap_data[$i]}")
    genome_size_aa["${viwrap_row_array[0]}"]="${viwrap_row_array[1]}"
    genome_completeness_value_aa["${viwrap_row_array[0]}"]="${viwrap_row_array[8]}"
    # Temporarily revert IFS changes to parse confidence method/value
    IFS=$whitespace_ifs
    read -a viwrap_completeness_method_value_array < <(echo "${viwrap_row_array[9]}")
    genome_completeness_confidence_method_aa["${viwrap_row_array[0]}"]="${viwrap_completeness_method_value_array[0]}"
    genome_completeness_confidence_value_aa["${viwrap_row_array[0]}"]="${viwrap_completeness_method_value_array[1]}"
    IFS=$'\037'
done

# Iterate through each row of the cluster report file
cluster_id=1
for cluster_report_row in "${cluster_report_data[@]}"; do
    echo "Processing cluster ${cluster_id} of ${total_clusters}"
    
    # Translate and iterate through the list of genomes in the cluster
    cluster_report_row_array=(${cluster_report_row})
    cluster_genomes_unparsed="${cluster_report_row_array[3]}"
    read -a cluster_genomes_array < <(echo $cluster_genomes_unparsed | tr ',' '\037' )
    
    # Set the cluster's maximally complete genome to the first in the array
    cluster_max_genome_name="${cluster_genomes_array[0]}"
    
    # Track the conditions under which the cluster's representative genome is chosen
    cluster_max_condition="DEFAULT"
    cluster_max_condition_set_flag=0

    # Compare remaining genomes to determine the cluster's maximally complete genome
    for i in "${!cluster_genomes_array[@]}"; do
        # Skip processing the first element in the array, since it is already being used as the default
        if (( i == 0 )); then
            continue
        fi

        # Assign genome name and completeness values to variables
        genome_name="${cluster_genomes_array[$i]}"
        genome_completeness_value="${genome_completeness_value_aa[$genome_name]}"

        # Handle genomes with no completeness value
        # Skip further processing for this genome if it's completeness value is empty
        # If no genomes in the cluster have a completeness value so far, compare genomes based on genome size
        if [ "$genome_completeness_value" = "" ]; then
            if [ "${genome_completeness_value_aa[$cluster_max_genome_name]}" = "" ]; then
                if [[ "${genome_size_aa[$genome_name]}" -gt "${genome_size_aa[$cluster_max_genome_name]}" ]]; then
                    cluster_max_genome_name=$genome_name 
                fi
                if (( cluster_max_condition_set_flag == 0 )); then
                    cluster_max_condition="DEFAULT_SIZE"
                    cluster_max_condition_set_flag=1
                fi
            fi
            continue
        # Check if the current genome is the only one in the cluster found to have a completeness value so far, setting it as the cluster's maximally complete genome
        elif [ "${genome_completeness_value_aa[$cluster_max_genome_name]}" = "" ]; then
            cluster_max_genome_name=$genome_name
            cluster_max_condition="COMPLETENESS"
            cluster_max_condition_set_flag=1
            continue
        fi
        
        # Compare the genome's completion to the cluster's current max completion found so far
        if (( $(echo "${genome_completeness_value_aa[$genome_name]} > ${genome_completeness_value_aa[$cluster_max_genome_name]}" | bc -l) )); then
            cluster_max_genome_name=$genome_name
        # When completeness values are equal, compare confidence levels
        elif (( $(echo "${genome_completeness_value_aa[$genome_name]} == ${genome_completeness_value_aa[$cluster_max_genome_name]}" | bc -l) )); then    
            cluster_max_genome_completeness_confidence_string="${genome_completeness_confidence_value_aa[$cluster_max_genome_name]}"
            genome_completeness_confidence_string="${genome_completeness_confidence_value_aa[$genome_name]}"
            if [[ "${enum_confidence_array[$genome_completeness_confidence_string]}" -gt "${enum_confidence_array[$cluster_max_genome_completeness_confidence_string]}" ]]; then
                cluster_max_genome_name=$genome_name     
            # When confidence levels are equal, compare genome sizes
            elif [[ "$genome_completeness_confidence_value" -eq "$cluster_max_completeness_confidence" ]]; then
                genome_size="${genome_size_aa[$genome_name]}"
                if [[ "${genome_size_aa[$genome_name]}" -gt "${genome_size_aa[$cluster_max_genome_name]}" ]]; then
                    cluster_max_genome_name=$genome_name                 
                fi
                cluster_max_condition="SIZE"
                cluster_max_condition_set_flag=1
            fi
            if (( cluster_max_condition_set_flag == 0 )); then
                cluster_max_condition="CONFIDENCE"
                cluster_max_condition_set_flag=1
            fi
        fi
        if (( cluster_max_condition_set_flag == 0 )); then
            cluster_max_condition="COMPLETENESS"
            cluster_max_condition_set_flag=1
        fi
    done

    # Concatenate tabulated data line and append to output file
    cluster_rep_data=""
    cluster_rep_data+="${cluster_id}\t"
    cluster_rep_data+="${#cluster_genomes_array[@]}\t"
    cluster_rep_data+="${cluster_max_genome_name}\t"
    cluster_rep_data+="${genome_completeness_value_aa[$cluster_max_genome_name]}\t"
    cluster_rep_data+="${genome_completeness_confidence_value_aa[$cluster_max_genome_name]}\t"
    cluster_rep_data+="${genome_completeness_confidence_method_aa[$cluster_max_genome_name]}\t"
    cluster_rep_data+="${genome_size_aa[$cluster_max_genome_name]}\t"
    cluster_rep_data+=$cluster_max_condition
    echo -e "${cluster_rep_data}" >> $tmp_output
    
    # Increment cluster_id value for the next row of 'CLUSTER_REPORT.txt'
    ((cluster_id=cluster_id+1))
done
echo "All clusters processed"

# Create output
cp $tmp_output $3
```

### Determine novelty of vOTUs
To determine whether any of the 1,865 unique vOTUs obtained from ViWrap were novel (i.e., not previously reported), we compared them using `skani` to all 42,811 "Bacteriophage" genomes on NCBI Virus (as of April 2024) and five different viral databases: Oral Virome Database (OVD, 48,425 genomes; PMID:35663034), Gut Virome Database (GVD, 33,242 genomes; PMID: 32841606), Gut Phage Database (GPD, 142,809 genomes; PMID: 33606979), Metagenomic Gut Virus catalogue (189,680 genomes, PMID: 34168315), and IMG/VR4 (15,722,824 genomes, PMID: 36399502).

```sh
# Sketched both databases and our vOTUs
skani sketch -l querylist.txt -t 8 -o query-sketch

# Compare .sketch files
skani dist -t 8 --ql samples.txt --rl refs.txt -o output
```
478 of the vOTUs from this study did not have a match at 95% ANI and 85% AF in these databases, indicating that these vOTUs have likely not been previously described.


## Human DNA Virus Identification

We also examined the metagenomic assemblies for human DNA virus genomes from the [NCBI Virus portal](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). The NCBI viral genomes FASTA file (downloaded on March 20th, 2024) had to be split into separate files in order for the genomes to be dereplicated via `anvi-dereplicate-genomes`, which was done with the following script:

#### split_fasta_file.sh

```sh
#!/bin/bash

# Splits a file of genomes downloaded from NCBI's virus genome database into individual files,
# where the name of the file is the virus' identifier
# Arguments:
#   1: NCBI compiled genome .fasta
#   2: Output directory

total_genomes=$(grep ">" $1 | wc -l)

count=0
while read line; do
    if [[ ${line:0:1} == ">" ]]; then
        ((count=count+1))
        read -a fasta_line_array < <(echo $line)
        fasta_file_identifier="${fasta_line_array[0]}"
        fasta_file_name="${fasta_file_identifier:1}.fasta"
        echo "Generating file ${fasta_file_name} (Genome ${count}/${total_genomes})"
        fasta_file_path="${2}/${fasta_file_name}"
        touch $fasta_file_path       
    fi
    echo "$line" >> $fasta_file_path
done < $1

echo "All genomes split"
```

The human virus genomes were then dereplicated with `anvi-dereplicate-genomes` at 95% ANI and 85% AF, as previously shown in the [vMAG dereplication step](#dereplication). Unlike the vMAG genome clusters, the dereplicated human virus genome clusters were not reassessed for new cluster-representative genomes, and the Anvi'o default of "centrality" was used select to the vOTUs.

### Mapping using BWA-MEM

BWA-MEM was then used through the `anvio-dev` conda environment to map the original metagenomic reads to the human DNA virus vOTUs. The database FASTA file used as input must first be indexed by BWA-MEM:

```sh
bwa index cluster_rep_human_database_fixed.fasta
```

```sh
conda run -p "miniforge3/envs/anvio-dev" bwa mem \
        -t 24 \
        -o "${sample_id}_bwamem_output.sam" \
        "$cluster_rep_human_database_fixed.fasta" \
        "$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_1.fastq)" \
        "$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_2.fastq)"
```

The mapping read counts of the `.sam` files from the BWA-MEM output were then extracted using the following shell and perl scripts:

#### extract_sam_counts.sh

```sh
#!/bin/bash

# Extracts the mapping counts from the contents of a .sam file and places them in a .txt file.
# Arguments:
#   1: Path to directory containing .sam files
#   2: Path to the perl script ("bwa_reference_summary_frag.pl")

for file in $1/*.sam; do
    file_name_prefix=$(basename $file .sam | cut -d_ -f1)
    perl $2 -f $file > "${1}/${file_name_prefix}_bwamem_counts.txt"
done
```

#### bwa_reference_summary_frag.pl

*Written by R. A. Richter at the J. Craig Venter Institute. (Present affiliation is the UC San Diego School of Medicine)*

```perl
#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Getopt::Long qw(GetOptionsFromArray);

my @files;
my $rate_flag;
my $verbose;
my $help;
my $status = GetOptionsFromArray(\@ARGV, 
				 "files=s{,}" => \@files,
				 #"sam" => \$sam_flag,
				 "rate" => \$rate_flag,
				 "verbose" => \$verbose,
				 "help" => \$help);

my $usage = 
    "This program takes BWA mapping files (SAM/BAM) and counts the number of mapped reads for each reference sequence.\n\n" .
    "Usage:$0\n" .
    "\t-f,--file\trequired\tstring(s)\tmapping-file1,...\n" .
    "\t-r,--rate\toptional\tflag\tratio if mulitple best hit?\t(default:off)\n" .
    "\t-v,--verbose\toptional\tflag\tverbose\t(default:off)\n" .
    "\t-h,--help\toptional\tflag\tthis message\n" .
    " e.g.:$0 -f test.bam -r\n";

if ( $help || @files == 0 ) {
    print $usage;
    exit;
}


my %Counts;
my %Sbjcts;		    ## subject matches for each query sequence
my $bwa_mem = 0;
my $first_read = 1;
foreach my $file ( @files ) {
    print STDERR "Loading $file\n";
    my $fh;
    ## Ignore unmapped reads (0x4) and supplementray alignments (0x800)
    my $exe = "samtools view -h -F 0x4 $file | samtools view -h -F 0x800 - |";
    open($fh, $exe) || die("Open error:$file");

    my $line = 0;
    my ( $prev, $curr );
    while ( <$fh> ) {
	## Get reference ids
	if ( /^\@SQ\s+SN:([^\s]+)/ ) {
	    $Counts{$1} = 0; ## Initialize all references match count to be 0
	} elsif ( /\@PG/ ) {
	    $bwa_mem = 1 if /mem/;
	}
	next if /^\@/;

	$line++;		
	parseLine( $_ );
    }
   
    print STDERR "#line:$line\n" ;
}

print STDERR "#read:" . scalar keys %Sbjcts, "\n";

## For each query sequence, increment match count to reference sequences or 
## evenly distribute to reference sequences (1/N).
foreach my $q ( keys %Sbjcts  ) {
    my @sbj = @{$Sbjcts{$q}};
    foreach my $s ( @sbj ) {
	if ( ! defined $rate_flag ) {
	    $Counts{$s}++;
	} else {
	    $Counts{$s} += 1/scalar(@sbj);
	}
    }
	
} 

foreach my $r ( sort keys %Counts ) {
    print join("\t", $r, $Counts{$r}), "\n";
}

sub getScore
{
    my ( $column ) = @_;
    $column =~ /NM:i:(\d+)/;
    return $1;
}

sub getOtherSbjcts
{
    my ( $edit, $Col) = @_;
    my @Sbj;
    ##AFNN01000024,-14652,101M,0
    if ( $Col->[$#$Col] =~ /XA:Z:(.+);$/ ) {
	my @others = split /;/, $1;
	foreach my $o ( @others ) {
	    my @f = split /,/, $o;
	    next if $f[$#f] > $edit;
	    push @Sbj, $f[0];
	}
    }
    return @Sbj;
}

sub parseLine
{
    my ($line) = @_;
    my @c = split /\s+/, $line;
	
    my $qid = $c[0];
    my $ref = $c[2];
    my $ecol = $bwa_mem ? 11 : 12;
    if ( $first_read ) {
	checkFormat($c[$ecol]);
	$first_read = 0;
    }
    my $edit = getScore($c[$ecol]);

    push @{$Sbjcts{$qid}}, $ref;
    push @{$Sbjcts{$qid}}, getOtherSbjcts( $edit, \@c );
}

sub checkFormat
{
    my ( $column ) = @_;
    $column =~ /NM:i:(\d+)/;
    die "Edit distance column error ($column)" unless $column =~ /^NM/;
}
```

This yields the amount of reads per node in the metagenomic assembly as a text file, and so the total reads associated with a given genome needed to be consolidated:

```sh
# Remove unique, contig-specific portions of contig deflines
for datafile in *_bwamem_counts.txt; do
    sed -e 's/_NODE.*\t/\t/' "$datafile" > "$datafile".fixed
done

# Add up all lines with the same field #1 (i.e. merge the counts of all contigs from each genome)
for datafile in *.fixed; do
    awk '$1!=p{ if (NR>1) print p, s; p=$1; s=0} {s+=$2} END{print p, s}' "$datafile" > "$datafile".totalled
done
```

The files containing the counts were then merged together to create the human DNA virus OTU table:

```sh
# Make OTU table
awk 'BEGIN {print "genome"} {print $1}' 1.count > 1st.txt

for datafile in *.totalled; do
    awk 'FNR ==1 {print FILENAME} {print $2}' "$datafile" > "$datafile.column"
done

paste 1st.txt *.column > otu_table.txt
```

Human DNA virus features with a sum of less than 10,000 reads were filtered out of the table.

## Final OTU Table Construction

We mapped our original metagenomic reads to both the phage and human DNA virus vOTUs to produce an OTU table for diversity analysis.

### Mapping with CoverM

[CoverM](https://github.com/wwood/CoverM) v0.6.1, implementing minimap2, was used to contruct an OTU table with relative abundances normalized for genome length:

```sh
conda run -p "ViWrap_conda_environments/ViWrap-Mapping" coverm genome \
        --threads 24 \
        --output-file "/sample_data/${sample_id}_coverm_output.txt" \
        --genome-fasta-directory "/all_cluster_rep_genomes/" \
        --genome-fasta-extension fasta \
        -1 "$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_1.fastq)" \
        -2 "$(find virome/${sample_id}_S*_R1_001_kneaddata_paired_2.fastq)"
```

CoverM created a relative abundance table for each sample. These files were joined together to create a relative abundance OTU table using the following script:

#### coverm_to_otu.sh

```sh
#!/bin/bash

# Compiles the output files of 'coverm genome' from a single directory into one OTU table file.
# Keeps relative abundance values and unmapped percentage.
# Arguments:
#   1: The directory of 'coverm genome' output files
#   2: Output file path, including file name

# Get file name
output_file_name=$(basename $2) 

# Set up temporary files
tmp_dir=$(mktemp -d)
tmp_coverm_column=$(mktemp -p $tmp_dir)
tmp_output=$(mktemp -p $tmp_dir)
tmp_paste=$(mktemp -p $tmp_dir)
trap 'rm -rf -- "$tmp_dir"' EXIT

# Build the OTU table from CoverM data
first_file_flag=0
for file in $1/*_coverm_output.txt; do
    if (( first_file_flag == 0 )); then
        head -n +1 $file > $tmp_output
        tail -n +2 $file | sort --stable -k1,1 >> $tmp_output
        first_file_flag=1
    else
        head -n +1 <(awk '{print $2}' $file) > $tmp_coverm_column
        awk '{print $2}' <(sort --stable -k1,1 <(tail -n +2 $file)) >> $tmp_coverm_column
        paste -d $'\t' $tmp_output $tmp_coverm_column > $tmp_paste
        cp $tmp_paste $tmp_output
    fi
done

cp $tmp_output $2
```

The values in the table were multiplied by a factor of 1,000,000 and rounded to the nearest integer using Microsoft Excel in order to process the OTU table as a frequency table for downstream analyses.

### Filtering

The final OTU table was filtered using `qiime feature-table filter-features`, to remove features with less than 100,000 total read counts across all samples (i.e., 0.1% total relative abundance) or that were present in less than 10 samples:

```sh
qiime feature-table filter-features --i-table otu_table.qza --p-min-frequency 100000 --p-min-samples 10 --o-filtered-table qiime_filtered_otu_table_1.qza
```

These filtering thresholds left us with 476 total features, comprising 475 phage vOTUs and 1 human DNA virus vOTU, *Human betaherpesvirus 7*.

### Looking for *Human betaherpesvirus 7* and *Human gammaherpesvirus 4* features
The MetaPhlAn2 analysis in the previous study (PMID:33239396) had identified the presence of *Human betaherpesvirus 7* (HHV7) and *Human gammaherpesvirus 4* (HHV4), and had suggested that both taxa were associated with caries, based on Songbird, with HHV4 being the taxa most associated with caries.  In contrast, HHV4 was not detected at significant abundances in the analysis performed in this study. Meanwhile, HHV7 was the only human DNA virus with a high enough relative abundance to make it into the final OTU table.  5_vRhyme_unbinned_20, obtained from ViWrap, was closely related to HHV7, but was a small ~5000bp fragment only. To more closely examine any HHV sequences in our dataset, the reads in Sample 5 (the sample with the highest abundance of both HHV7 and HHV4 in the previous study) were mapped to the annotated HHV7 and HHV4 genomes:

```sh
ncbi-genome-download viral -F gff -A GCF_002402265.1

ncbi-genome-download viral -F fasta -A GCF_002402265.1

minimap2 -ax sr GCF_002402265.1_ASM240226v1_genomic.fna ../5_S4_R1_001_kneaddata_paired_1.fastq ../5_S4_R1_001_kneaddata_paired_2.fastq > hhv4_sample_5.sam

samtools sort -o hhv4_sample_5.sorted.bam hhv4_sample_5.sam

samtools index hhv4_sample_5.sorted.bam

featureCounts -a GCF_002402265.1_ASM240226v1_genomic.gff -o counts -t gene -g ID -p  hhv4_sample_5.sorted.bam
```
Only 72 reads mapped to the HHV4 genome from Sample 5 using this method. In contrast, 8,112 reads mapped to the HHV7 genome from Sample 5, with reads mapping generally evenly throughout the genome. The ViWrap vOTU with similarity to HHV7, 5_vRhyme_unbinned_20, mapped to the 4210-9272 and 147295-152318 (terminal repeat regions) of RefSeq HHV7 with 99% identity and did NOT match HHV4. We also re-examined Sample 5 using an updated version of MetaPhlAn, v4.1.0, with the `--mpa3` and `--add_viruses` flags:
```sh
metaphlan 5_S4_R1_001_kneaddata_paired_1.fastq.bowtie2out.txt --input_type bowtie2out --nproc 12 --mpa3 --add_viruses > 5_metaphlan_profile_virus.txt
```
 This did not detect a significant abundance of either HHV7 or HHV4 using the default cutoffs (as was employed in the previous study). It is unclear why the previous study detected significant levels of HHV4 (mainly in Sample 5), which were not detected by any methods here. The most likely explanation would seem to be differences in the MetaPhlAn algorithm and/or marker gene databases between the current version and those from ~2018, when the previous analyses were performed.


## Downstream Data Analysis

With the final OTU table constructed and filtered, downstream analyses of the caries-associated virome were performed.

### QIIME 2 Diversity Tests

We utilized QIIME 2's inbuilt diversity analysis tool, `qiime diversity`, to test different diversity metrics, and any significant association with sample metadata.

#### Relative Frequency Table

```sh
# Relative abundance table for qiime pcoa-biplot
conda run -p $2 qiime feature-table relative-frequency \
    --i-table $1 \
    --o-relative-frequency-table $file_otu_relative_frequency
```

#### Alpha Diversity

```sh
# Alpha diversity
conda run -p $2 qiime diversity alpha \
    --i-table $1 \
    --p-metric shannon \
    --o-alpha-diversity $file_alpha_diversity

# Alpha correlation spearman
conda run -p $2 qiime diversity alpha-correlation \
    --i-alpha-diversity $file_alpha_diversity \
    --m-metadata-file $file_sample_metadata \
    --p-method spearman \
    --o-visualization $file_alpha_correlation_spearman
```

#### Beta Diversity

```sh
# Beta diversity braycurtis distance matrix
conda run -p $2 qiime diversity beta \
    --i-table $1 \
    --p-metric braycurtis \
    --o-distance-matrix $file_beta_diversity_braycurtis_distance_matrix

# Beta diversity braycurits pcoa
conda run -p $2 qiime diversity pcoa \
    --i-distance-matrix $file_beta_diversity_braycurtis_distance_matrix \
    --o-pcoa $file_beta_diversity_braycurtis_pcoa

# Beta diversity braycurtis biplot
conda run -p $2 qiime diversity pcoa-biplot \
    --i-pcoa $file_beta_diversity_braycurtis_pcoa \
    --i-features $file_otu_relative_frequency \
    --o-biplot $file_beta_diversity_braycurtis_biplot

# Beta diversity braycurtis permanova
conda run -p $2 qiime diversity beta-group-significance \
    --i-distance-matrix $file_beta_diversity_braycurtis_distance_matrix \
    --m-metadata-file $file_sample_metadata \
    --m-metadata-column Status \
    --p-method permanova \
    --o-visualization $file_beta_diversity_braycurtis_permanova
```
*The above steps were also performed using the Jaccard metric for beta diversity.*

We used the QIIME 2 plugin [DEICODE](https://library.qiime2.org/plugins/deicode/19/) to create a robust Aitchison PCA (RPCA) distance matrix:

```sh
# Beta diversity rpca
conda run -p $2 qiime deicode rpca \
    --i-table $1 \
    --p-min-feature-count 10 \
    --p-min-sample-count 500 \
    --o-biplot $file_beta_diversity_rpca_biplot \
    --o-distance-matrix $file_beta_diversity_rpca_distance_matrix

# Beta diversity rpca permanova
conda run -p $2 qiime diversity beta-group-significance \
    --i-distance-matrix $file_beta_diversity_rpca_distance_matrix \
    --m-metadata-file $file_sample_metadata \
    --m-metadata-column Status \
    --p-method permanova \
    --o-visualization $file_beta_diversity_rpca_permanova
```

The visualization of these biplots was done using the QIIME 2 [Emperor](https://biocore.github.io/emperor/) plugin:

```sh
# braycurtis
conda run -p $2 qiime emperor biplot \
    --i-biplot $file_beta_diversity_braycurtis_biplot \
    --m-sample-metadata-file $file_sample_metadata \
    --m-feature-metadata-file $file_feature_metadata \
    --o-visualization $file_beta_diversity_braycurtis_biplot_visualization \
    --p-number-of-features 15

# jaccard
conda run -p $2 qiime emperor biplot \
    --i-biplot $file_beta_diversity_jaccard_biplot \
    --m-sample-metadata-file $file_sample_metadata \
    --m-feature-metadata-file $file_feature_metadata \
    --o-visualization $file_beta_diversity_jaccard_biplot_visualization \
    --p-number-of-features 15

# rpca
conda run -p $2 qiime emperor biplot \
    --i-biplot $file_beta_diversity_rpca_biplot \
    --m-sample-metadata-file $file_sample_metadata \
    --m-feature-metadata-file $file_feature_metadata \
    --o-visualization $file_beta_diversity_rpca_biplot_visualization \
    --p-number-of-features 15
```

An annotated version of the RPCA biplot is featured as Figure 2A.

### Differential Abundance
Due to the issues inherent with performing differential abundance analysis on compositional data (Gloor et al., 2017, PMID:29187837; Morton et al., 2019, PMID:31222023), we used three different methods to examine differential abundance: Songbird, DESeq2, and ANCOM-BC.

#### 1. Songbird

The QIIME 2 plugin [Songbird](https://github.com/biocore/songbird) v1.0.4 was used to generate a table of each feature's differential ranking.

```sh
qiime songbird multinomial \
    --i-table qiime_filtered_otu_table_1.qza \
    --m-metadata-file metadata.txt \
    --p-formula "Status" \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --p-epochs 100000 \
    --verbose \
    --o-differentials differentials.qza \
    --o-regression-stats regression_stats.qza \
    --o-regression-biplot regression_biplot.qza
```

Using Songbird for a dataset also necessitates the generation of a [null model](https://github.com/biocore/songbird?tab=readme-ov-file#612-null-models-and-qiime-2--songbird).

```sh
# Null model (identical input except for formula and output files)
qiime songbird multinomial \
    --i-table qiime_filtered_otu_table_1.qza \
    --m-metadata-file metadata.txt \
    --p-formula "1" \
    --p-differential-prior 0.5 \
    --p-summary-interval 1 \
    --p-epochs 100000 \
    --verbose \
    --o-differentials null_differentials.qza \
    --o-regression-stats null_stats.qza \
    --o-regression-biplot null_biplot.qza

# Visualize the original and null models together
qiime songbird summarize-paired \
    --i-regression-stats regression_stats.qza \
    --i-baseline-stats null_stats.qza \
    --o-visualization paired_summary.qzv
```

The Songbird output file `differentials.qza` was converted into the metadata file `virus_genome_metadata.txt` for our viral genomes:

```sh
qiime tools export --input-path differentials.qza --output-path ./

head -n +1 differentials.tsv > virus_genome_metadata.txt && tail -n +3 differentials.tsv | sort -k 3n | column -s $'\t' -t >> virus_genome_metadata.txt
```

#### 2. DESeq2
The following R code was used to determine the vOTUs with significantly different relative abundances between the caries and healthy groups:

```sh
## DA analysis with DESeq2 
## Updated for Oral Virome 4/10/24

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("genefilter")

# Set working directory
setwd("~/Desktop/Virome/")
getwd()

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts

countdata <- read.table("qft_1_otu_table.txt", header=TRUE, row.names=1)

# Convert to matrix
#countdata = as.matrix(UAB_deSeq_input)
countdata <- as.matrix(countdata)
head(countdata)

# Assign samples
(condition <- factor(c(rep("caries", 22), rep("healthy", 23))))

# DESeq analysis----------------------------------------------------

library(DESeq2)


# Create a coldata frame and instantiate the DESeqDataSet 
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

#set healthy as reference level
dds$condition <- relevel(dds$condition, "healthy")

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots 
## Use RColorBrewer
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(11, 11), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition") + theme_bw()
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-50, 50))
dev.off()


#  DE results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized counts
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Viral_genome"
head(resdata)

## Write results
write.csv(resdata, file="deseq2-results.csv")

## Examine p-values
hist(res$pvalue, breaks=50, col="grey")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Viral_genome, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Viral_genome, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()


#Heat map for the most abundant vOTUss

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds
dds <- estimateSizeFactors(dds)
normalizeddata <- counts(dds, normalized=TRUE)
write.csv(normalizeddata, file="normalizeddata.csv")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition")]) 
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE)
```
#### 3. ANCOM-BC

ANCOM-BC was used through QIIME 2, and the resulting differentials were exported and visualized.

```sh
qiime composition ancombc \
    --i-table qiime_filtered_otu_table_1.qza \
    --m-metadata-file metadata.tsv \
    --p-formula Status \
    --o-differentials ancombc_differentials.qza
```

### MMvec

Using feature tables from the previous study, the QIIME 2 plugin [MMvec](https://github.com/biocore/mmvec/tree/master) v1.0.6 was used to obtain conditional probabilities of co-occurrences in two different multiomics datasets: the viral genomes and human host cytokine levels, and the viral genomes and bacterial levels.

The bacterial OTU table was from [Supplemental Table S3 of the previous study](https://genome.cshlp.org/content/suppl/2020/12/16/gr.265645.120.DC1/Supplemental_Table_S3.xlsx), and was converted to the QIIME 2-compatible file `filtered_bacteria_otu_table.qza` with the following commands after being exported to a tab-delimited format in Microsoft Excel:

```sh
# Activate conda environment
conda activate qiime2-2020.6

# Exclude the table descriptor header
tail -n +2 Baker_Supplemental_Table_S3.txt > bacteria_otu_table.txt

# Remove the trailing space of the species names
sed 's/ //g' bacteria_otu_table.txt
```

The OTU table was then filtered, such that it only contains features that also have metadata in the corresponding metadata file [Supplemental Table S4](https://genome.cshlp.org/content/suppl/2020/12/16/gr.265645.120.DC1/Supplemental_Table_S4.xlsx) of the previous study. This was done with the following Python script to generate `filtered_bacteria_otu_table.txt`:

#### filter_bacteria_frequency_table.py

```py
import csv

def import_tsv_data(file_name: str, data_structure: str = "list", header: bool = True) -> tuple:
    if not header:
        file_fields_list = None

    if data_structure == "list":
        file_data = []
    elif data_structure == "dict":
        file_data = {}

    fields_line_flag = True
    with open(file_name, 'r') as file:
        for line in file:
            line_list = line.split('\t')
            line_list[-1] = line_list[-1].strip()
            if fields_line_flag and header:
                file_fields_list = line_list
                fields_line_flag = False
            else:
                if data_structure == "list":
                    file_data.append(line_list)
                elif data_structure == "dict":
                    file_data[line_list[0]] = line_list

    return file_fields_list, file_data


bacteria_otu_fields, bacteria_otu_data = import_tsv_data("bacteria_otu_table.txt", "list")
bacteria_metadata_fields, bacteria_metadata_data = import_tsv_data("bacteria_metadata.txt", "dict")

print(bacteria_metadata_data)

with open("filtered_bacteria_otu_table.txt", 'w', newline='') as file_output:
    csv_writer = csv.writer(file_output, delimiter='\t')
    csv_writer.writerow(bacteria_otu_fields)
    for frequency_list in bacteria_otu_data:
        if frequency_list[0] in bacteria_metadata_data:
            csv_writer.writerow(frequency_list)
```

The filtered table was then imported into QIIME 2:

```sh
# Convert .txt to .biom file
biom convert -i filtered_bacteria_otu_table.txt -o filtered_bacteria_otu_table.biom --to-hdf5

# Import .biom file to QIIME 2 as a Frequency FeatureTable artifact
qiime tools import --input-path filtered_bacteria_otu_table.biom --output-path filtered_bacteria_otu_table.qza --type FeatureTable[Frequency]
```

The `bacteria_metadata.txt` metadata file was converted from the previous study's Supplemental Table S4 with the following commands after being exported to a tab-delimited format in Microsoft Excel:

```sh
tail -n +2 Baker_Supplemental_Table_S4.txt > bacteria_metadata.txt
```

The `cytokines.txt` frequency table and `cytokines_metadata_ACTUAL.txt` files were imported to QIIME 2 in a pipeline analogous to the above, with the file names `cytokines.qza` and `cytokines_metadata_ACTUAL.txt`, respectively.

In both the bacteria and cytokine metadata files, `vim` was used to edit the header of the features to '`id`'.

The virus and bacteria frequency tables were input to MMvec to generate the conditional probabilities of co-occurrences:

```sh
qiime mmvec paired-omics --i-microbes qiime_filtered_otu_table_1.qza --i-metabolites /bacteria_otu_table.qza --p-summary-interval 1 --output-dir mmvec_output_1
```

The same was done for the conditial probabilities of co-occurrences of the viruses and cytokines; to generate the biplot where the viruses are represented as the points and the cytokines as the vectors, the `--i-microbes` and `--i-metabolites` arguments were switched, like so:

```sh
qiime mmvec paired-omics --i-metabolites /qiime_filtered_otu_table_1.qza --i-microbes cytokines.qza --p-summary-interval 1 --p-input-prior 0.1 --output-dir mmvec_output_1
```