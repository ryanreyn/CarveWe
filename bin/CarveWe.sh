#!/bin/bash
#This file is the main file that we will build out to incorporate submodules of our
#genome annotation, model and media prediction, and growth sensitivity pipeline

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
    printf "  a variety of the submodules developed for the CarveWe pipeline. Currently this\n"
    printf "  version supports the alignment of input genomes to the genomes used to build\n"
    printf "  the SOM neural network described in the associated publication to this work\n"
    printf "  via best hit BLAST alignment. Additional functionality to come!\n\n"

    printf "\n ---------------------------- ${LIGHTYELLOW}REQUIRED INPUTS${NC} ----------------------------\n\n"
    printf "    1) Input the path to a directory of FASTA formatted files or to an individual file\n"
    printf "        - [-d <directory>] default: false\n"
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
    printf "    ${LIGHTYELLOW}Number of threads for BLAST:${NC}\n\n"
    printf "        - [-t <int>] default: 1\n"
    printf "            Specify the number of threads for BLAST multi-threading to use\n\n"

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

#Setting up and parsing arguments
#Setting some default parameters
out_dir='CarveWe_output'
fasta_dir='false'
num_threads=1
fna_type='false'
faa_type='false'

#Establish the options for this program
while getopts p:dano:t: args
do
    case ${args} in
        p) fasta_file=${OPTARG};;
        d) fasta_dir='true';;
        a) faa_type='true';;
        n) fna_type='true';;
        o) out_dir=${OPTARG};;
        t) num_threads=${OPTARG};;
        \?) printf "\n  ${RED}Invalid argument: -${OPTARG}${NC}\n\n    Run 'CarveWe' with no arguments, '-h', or '--help' only to see help menu.\n\n" >&2 && exit
    esac
done

###### This section is checking and reporting which blastdb files will be used ######
blastdb_dir="../blastdb_files" #This needs to be adjusted when this is released as 
#a package.
printf " ${YELLOW}BLAST database files to be used:${NC}\n\n"

if [ "$fna_type" = 'true' ]
then
    blastdb_file="$blastdb_dir/carvewe-hq-genomes"
    printf " ${YELLOW}carvewe-hq-genomes. This is a database file of the nucleotide${NC}\n"
    printf " ${YELLOW}sequences from all of the CarveWe genomes${NC}\n\n"
elif [ "$faa_type" = 'true' ]
then
    blastdb_file="$blastdb_dir/carvewe-protein-seqs"
    printf " ${YELLOW}carvewe-protein-seqs. This is a database file of the amino acid${NC}\n"
    printf " ${YELLOW}sequences from all of the CarveWe genomes${NC}\n\n"
fi

#Call the first subprocess, conditioning on whether the -n or -a flag was provided
printf " ${YELLOW}Running the BLAST subprocess to align input genomes using best hit.${NC}\n\n"
if [ "$fna_type" == true ]
then
    genome-blast.sh -p ${fasta_file} -d ${fasta_dir} -n \
    -o ${out_dir} -t ${num_threads} -b ${blastdb_file}
elif [ "$faa_type" == true ]
then
    genome-blast.sh -p ${fasta_file} -d ${fasta_dir} -a \
    -o ${out_dir} -t ${num_threads} -b ${blastdb_file}
fi

