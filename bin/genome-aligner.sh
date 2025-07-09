#!/bin/bash
#This is a process that is also integrated into the CarveWe tool enables a user to
#provide their own genome(s) and get the closest hit to a SOM cluster via BLAST
#Establish some flags for reading user input and providing help

#Setting some colors to use for printing verbose output (inspired by Mike Lee)
GREEN='\e[32m'
RED='\e[31m'
YELLOW='\e[33m'
LIGHTYELLOW='\e[93m'
NC='\e[0m'

############## Defining help info ###################################################
help_menu()
{
    printf "\n ---------------------------- ${LIGHTYELLOW}HELP INFO${NC} ---------------------------- \n\n"
    printf "  This program accepts input genomes in the form of a path to a directory of\n"
    printf "  fasta files or the path to an individual fasta file of interest and performs\n"
    printf "  a best hit BLAST alignment to the genomes used to build the SOM neural network.\n"
    printf "  This process supports the input of both nucleotide and protein fasta sequences.\n"

    printf "\n ---------------------------- ${LIGHTYELLOW}REQUIRED INPUTS${NC} ----------------------------\n\n"
    printf "    1) Input the path to a directory of FASTA formatted files or to an individual file\n"
    printf "        - [-d ] default: false\n"
    printf "            Define whether the input path is a file or directory\n"
    printf "        - [-p <file>] file path to an individual FASTA file or directory of files\n"
    printf "    2) Define whether input genomes are in nucleotide (.fna) or protein (.faa) format\n"
    printf "        This choice is mutually exclusive, both are by default set to false.\n"
    printf "        - [-n ] Specify that the FASTA file(s) are in nucleotide format\n"
    printf "        - [-a ] Specify that the FASTA file(s) are in amino acid format\n"

    printf "\n ---------------------------- ${LIGHTYELLOW}OPTIONAL INPUTS${NC} ----------------------------\n\n"
    printf "    ${LIGHTYELLOW}Output directory specification:${NC}\n\n"
    printf "        - [-o <str>] default: CarveWe_output\n"
    printf "            Specify the desired output directory\n\n"
    printf "    ${LIGHTYELLOW}Number of threads for parallelizable subprocesses:${NC}\n\n"
    printf "        - [-t <int>] default: 1\n"
    printf "            Specify the number of threads for any multi-threading steps\n\n"

    exit
}

#Help menu to be displayed when user inputs the command with no inputs or with
#the -h or --help arguments
if [ "$#" == 0 ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then
    help_menu
fi

#Check for all required program dependencies
if ! command -v blastn > /dev/null || ! command -v blastp > /dev/null
then
    printf "\n ${RED}BLAST is an essential dependency but does not appear to be in your PATH${NC}\n"
    printf "\n Install BLAST or add your exisiting BLAST installation to your PATH.\n\n"
    printf "\nExiting for now.\n\n"
fi

fasta_file=""
is_dir="true"
faa_type="false"
fna_type="false"
blastdb_dir="../ref_data"

while getopts "p:d:ano:t:b:" args; do
    case "${args}" in
        p) fasta_file=${OPTARG};;
        d) is_dir=${OPTARG};;
        a) faa_type='true';;
        n) fna_type='true';;
        o) out_dir=${OPTARG};;
        t) num_threads=${OPTARG};;
    esac
done

#Make the output directory
if [ -d ${out_dir} ]
then
    rm -r ${out_dir}
    mkdir ${out_dir}
else
    mkdir ${out_dir}
fi

#Do some initial checks
combined_fasta="${out_dir}/tmp_search.fasta"
touch ${combined_fasta}
if [ "$is_dir" != true ]
then
    echo "Provided with a file."
    cat ${fasta_file} > ${combined_fasta}
elif [ "$is_dir" == true ]
then
    echo "Provided with a directory."
    fasta_list="${out_dir}/tmp_fasta-list.txt"
    ls ${fasta_file} > ${fasta_list}
    while read line
    do
        cat ${fasta_file}'/'${line} >> ${combined_fasta}
    done < ${fasta_list}
fi

`head -n 10 ${combined_fasta} | echo`

#Initially building it using my conda infrastructure but will remove so all that's
#required is to have BLAST installed
module purge
eval "$(conda shell.bash hook)"

conda activate NGStools

outfile="${out_dir}blast-to-carvewe.tsv"

#Run the appropriate blastn or blastp command depending on the n/a flag passed to
#the top level process call
if [ "$fna_type" == true ]
then
    echo "Running blastn"
    blastdb_file="$blastdb_dir/carvewe-hq-genomes"
    blastn -num_threads ${num_threads} \
    -query ${combined_fasta} \
    -db ${blastdb_file} \
    -outfmt 6 \
    -perc_identity 95 \
    -qcov_hsp_perc 100 \
    -num_alignments 1 \
    -out ${outfile} 2>>$outdir"log.txt"
elif [ "$faa_type" == true ]
then
    echo "Running blastp"
    blastdb_file="$blastdb_dir/carvewe-protein-seqs"
    blastp -num_threads ${num_threads} \
    -query ${combined_fasta} \
    -db ${blastdb_file} \
    -outfmt 6 \
    -num_alignments 1 \
    -out ${outfile} 2>>$outdir"log.txt"
fi

#Post-process the output BLAST table to extract query genome, CarveWe genome
genome_file=${out_dir}'/tmp-genome-hits-file.txt'
`cut -f 2 ${outfile} | sed 's/-scaffold_[0-9]*$//g' > ${genome_file}`

genome_alignments=${out_dir}'/tmp-genome-alignments.tsv'
`cut -f 1,2 ${outfile} | sed 's/-scaffold_[0-9]*$//g' > ${genome_alignments}`

#Pull in cluster assignment file and align with processed BLAST output
cluster_assignment="../final-cluster-assignments.tsv"

match_file=${out_dir}'/tmp-cluster-matches.txt'
`grep -f ${genome_file} ${cluster_assignment} > ${match_file}`

#Merge files to get a final file of the query genome and the matched genomes
final_file=${out_dir}'/cluster-matches.tsv'
paste ${genome_alignments} ${match_file} > ${final_file}

#Remove any temporary files generated during the BLAST process
rm -f ${out_dir}/tmp*

#Deactivate conda environment before releasing back to main process
conda deactivate
