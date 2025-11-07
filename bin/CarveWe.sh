#!/usr/bin/env bash
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
    printf "    ${LIGHTYELLOW}Run BLAST against genomes used to build SOM clusters:${NC}\n\n"
    printf "        - [-b ] default: false\n"
    printf "            Specify whether you wish to run the BLAST comparison or not\n\n"
    printf "            The BLAST alignment typically takes about an hour to run\n\n"
    printf "    ${LIGHTYELLOW}Model ensemble size for the carving subprocess:${NC}\n\n"
    printf "        - [-e <int>] default: 2\n"
    printf "            Specify how many independent models for carve to generate for each genome\n\n"
    printf "    ${LIGHTYELLOW}Model ensemble size for the carving subprocess:${NC}\n\n"
    printf "        - [-s <int>] default: false\n"
    printf "            Specify whether to skip carveme model predictions\n\n"

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
# Check for CarveWe environment
if [ -z "${CARVEWE_INSTALL_DIR:-}" ]; then
    printf "\n ${RED}CarveWe environment variables not set. Is your conda environment activated?${NC}\n\n"
    exit 1
fi

#Setting some default parameters
out_dir='CarveWe_output'
fasta_dir='false'
num_threads=1
fna_type='false'
faa_type='false'
ensemble_size=2
run_blast='false'
data_dir="${CARVEWE_DATA_DIR}"
scripts_dir="${CARVEWE_SCRIPTS_DIR}"
skip_carving='false'

# Ensure scripts are found in PATH
export PATH="${scripts_dir}:${PATH}"

#Establish the options for this program
while getopts p:dano:t:e:bs args
do
    case ${args} in
        p) fasta_file=${OPTARG};;
        d) fasta_dir='true';;
        a) faa_type='true';;
        n) fna_type='true';;
        o) out_dir=${OPTARG};;
        t) num_threads=${OPTARG};;
        e) ensemble_size=${OPTARG};;
        b) run_blast='true';;
        s) skip_carving='true';;
        \?) printf "\n  ${RED}Invalid argument: -${OPTARG}${NC}\n\n    Run 'CarveWe' with no arguments, '-h', or '--help' only to see help menu.\n\n" >&2 && exit
    esac
done

#Normalizing directories to remove any lagging "/" symbols
out_dir=`echo "${out_dir%/}"`
data_dir=`echo "${data_dir%/}"`


#Initialize the log file for error output
touch $out_dir"/log.txt"

##### This section will only run if the user elects to run the BLAST subprocess ######
##### It will pull the necessary genome files from the web and build BLAST db ########
#Contain the process within a file check to see if the user has already generated the
#blast files
if [ "$run_blast" == "true" ]
then
    printf "${YELLOW}You have elected to run the BLAST, checking if the database files exist.${NC}\n\n"
    if [ "$fna_type" == "true" ] && [ ! -f  "$data_dir/carvewe-hq-genomes.nhr" ]
    then 
        printf "${RED}BLAST database files not found, so we will generate them.${NC}\n\n"
        "${scripts_dir}/construct-db.sh" -n -o $out_dir -p $num_threads
    elif [ "$faa_type" == "true" ] && [ ! -f "$data_dir/carvewe-protein-seqs.pdb" ]
    then
        printf "${RED}BLAST database files not found, so we will generate them.${NC}\n\n"
        "${scripts_dir}/construct-db.sh" -a -o $out_dir -p $num_threads
    fi
else
    printf "${YELLOW}You have already generated the BLAST database files!${NC}\n\n"
fi
###### This section will run the first subprocess to BLAST input genomes to SOM ######
###### genomes ############################################
#This first subprocess is optional so we will enclose it in a conditional based on
#user flag inputs
if [ "$run_blast" == "true" ]
then
    printf " ${YELLOW}BLAST database files to be used:${NC}\n\n"

    if [ "$fna_type" = 'true' ]
    then
        printf " ${YELLOW}carvewe-hq-genomes. This is a database file of the nucleotide${NC}\n"
        printf " ${YELLOW}sequences from all of the CarveWe genomes${NC}\n\n"
    elif [ "$faa_type" = 'true' ]
    then
        printf " ${YELLOW}carvewe-protein-seqs. This is a database file of the amino acid${NC}\n"
        printf " ${YELLOW}sequences from all of the CarveWe genomes${NC}\n\n"
    fi

    #Call the BLAST program, conditioning on whether the -n or -a flag was provided
    printf " ${YELLOW}Running the BLAST subprocess to align input genomes using best hit.${NC}\n\n"
    
    if [ "$fna_type" == true ]
    then
        "${scripts_dir}/genome-aligner.sh" -p ${fasta_file} -d ${fasta_dir} -n \
        -o ${out_dir} -t ${num_threads}
    elif [ "$faa_type" == true ]
    then
        "${scripts_dir}/genome-aligner.sh" -p ${fasta_file} -d ${fasta_dir} -a \
        -o ${out_dir} -t ${num_threads}
    fi
else
    printf "${YELLOW}Skipping the BLAST alignment step!${NC}\n\n"
fi

##### This section will run segments of the CarveWe pipeline for users interested in
#processing their own genomes through CarveMe, COBRApy, and metabolic sensitivity####
#Run the model carving and reaction info extraction subprocess
if [ "$skip_carving" == 'true' ]
then
    printf "${YELLOW}Skipping the model prediction step!${NC}\n\n"
else
    printf "${YELLOW}Running the provided genomes through CarveMe to annotate metabolic models:${NC}\n\n"

    if [ "$fna_type" == 'true' ]
    then
        "${scripts_dir}/run-carveme.sh" -n -i ${fasta_file} -o ${out_dir} -e ${ensemble_size} \
        -p ${num_threads}
    else
        "${scripts_dir}/run-carveme.sh" -a -i ${fasta_file} -o ${out_dir} -e ${ensemble_size} \
        -p ${num_threads}
    fi

    #Building a subprocess here to extract model quality and generate some analytical plots
    #for users to examine 
    xml_dir=$out_dir"/xml_files"
    if [[ "$ensemble_size" -ge 10 ]]
    then
        printf "${YELLOW}Passing model files to an ensemble quality prediction subprocess${NC}\n\n"
        python "${scripts_dir}/extract-model-quality.py" -d $xml_dir -w $out_dir
    else
        printf "${RED}The specified ensemble size is too small to reliably analyze within ensemble model differences${NC}\n\n"
    fi
fi

#Now we will pass the xml and reaction info files to the COBRApy subprocess
media_dir=$out_dir"/media_files"
if [ -d "$media_dir" ]
then
    rm -r $media_dir
    mkdir $media_dir
else
    mkdir $media_dir
fi

#Create a temporary file with all of the genome files (dropping the fasta extension)
printf "${YELLOW}\nRunning the SBML CarveMe model files through COBRApy for media prediction:${NC}\n\n"

genome_list=$out_dir"/tmp-genomes-list.txt"
ls $out_dir/xml_files/*.xml | sed "s/\.xml//g; s/.*\///g" > $genome_list

while read genome
do
    printf "Running $genome through media prediction\n\n"
    python "${scripts_dir}/generate_ensemble_media_microbiomics.py" --ensemble-size $ensemble_size \
    --work-dir $out_dir --media-dir $media_dir $genome

    printf "Completed running $genome through media prediction, outputting the results to $media_dir\n\n"
done < $genome_list

#Now we will pass the raw COBRApy predictions to a script to filter, merge, and 
#convert the output in order to run metabolite sensitivity tests as well
printf "${YELLOW}Converting raw COBRApy output to correct input for sensitivity prediction:${NC}\n\n"

python "${scripts_dir}/convert-media-output.py" --media-dir $media_dir --work-dir $out_dir


#Now we will pass our reformatted and combined media data to predict metabolite
#sensitivities for all input genomes
printf "${YELLOW}Extracting sensitivity values for our defined metabolite classes:${NC}\n\n"

python "${scripts_dir}/get_met_depends.py" --media-dir $media_dir --work-dir $out_dir \
    --ensemble-size $ensemble_size --data-dir $data_dir

#Remove temporary file listing out the genomes
rm $genome_list
