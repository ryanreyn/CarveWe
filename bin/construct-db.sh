#!/bin/bash
#This script will pull in the fasta files for the genomes used to build the SOM clusters
#We will set it up to grab these files from the remote repository, build a blast
#database, and then use this in the blast alignment. We will set this to delete after
#by default
fna_type="false"
faa_type="false"
delete_tmp="true"
max_processes=1

while getopts "nado:p:" args; do
    case "${args}" in
        n) fna_type='true';;
        a) faa_type='true';;
        d) delete_tmp='false';;
        o) out_dir=${OPTARG};;
        p) max_processes=$((OPTARG - 1));;
    esac
done

#Set up some initial variables
files="hq-genome-ids.txt"
tmp_dir=$out_dir"tmp/"
fna_blastdb="${CARVEWE_DATA_DIR}/carvewe-hq-genomes"
faa_blastdb="${CARVEWE_DATA_DIR}/carvewe-protein-seqs"

if [[ -d ${tmp_dir} ]]
then
    printf "Making directory!\n\n"
    rm -r ${tmp_dir}
    mkdir ${tmp_dir}
else
    mkdir ${tmp_dir}
fi

`cp ${CARVEWE_DATA_DIR}/${files} ${tmp_dir}`

here=`pwd`
cd ${tmp_dir}

#Create the initial piece of the web address to index
prefix="https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/1.0/genomes/"

if [ "$fna_type" == true ]
then
    concat_file="concat-fna-files.fna"
    touch $concat_file
elif [ "$faa_type" == true ]
then
    concat_file="concat-faa-files.faa"
    touch $concat_file
fi

#We will run the fasta file acquisitions using background multiprocessing to speed up
curr_processes=0

#We also want a progress bar rather than dumping all the standard output to terminal
curr_iter=1
max_iter=`wc -l $files | cut -d ' ' -f 1`

#Run the loop that will grab all of the necessary fasta files
while read genome
do
    while [[ $current_processes -ge $max_processes ]]
    do
        wait -n # Waiting for any process to complete to add more processes
        current_processes=$((current_processes - 1))
    done
    if [ "$faa_type" == true ]
    then
        (
        suffix=${genome}"/gene_call/"${genome}"-prokka.faa"
        wget --no-check-certificate -a "${out_dir}/log.txt" ${prefix}${suffix}
        ) &
        current_processes=$((current_processes + 1))
    elif [ "$fna_type" == true ]
    then
        (
        suffix=${genome}"/prokka/"${genome}"-Bacteria-prokka/"${genome}"-prokka-original.fna"
        wget --no-check-certificate -a "${out_dir}/log.txt" ${prefix}${suffix}
        ) &
        current_processes=$((current_processes + 1))
    fi
    #Adding a progress bar instead of letting standard wget output hit terminal
    # Calculate progress percentage
    progress_percent=$(( (curr_iter * 100) / max_iter ))

    # Build the progress bar
    bar_length=50
    completed_length=$(( progress_percent * bar_length / 100 ))
    remaining_length=$(( bar_length - completed_length ))
    completed_bar=$(printf "%${completed_length}s" | tr " " "#")
    remaining_bar=$(printf "%${remaining_length}s" | tr " " " ")

    # Print the progress bar
    printf "\rFasta file download progress: [%s%s] %3d%%" "$completed_bar" "$remaining_bar" "$progress_percent"

    # Increment the counter
    curr_iter=$((curr_iter + 1))
done < $files

#Make sure to wait for any processes that are still running to finish
while [[ $current_processes -gt 0 ]]
do
    wait -n
    current_processes=$((current_processes - 1))
done

printf "\nAll files downloaded!\n\n"

#Now we will concatenate all of the files together into one file for making the db
list_files="tmp.txt"
ls ./[A-Z]*.f[an]a > $list_files

curr_iter=1

printf "Concatenating fasta files!\n\n"
while read line
do
    cat $line >> $concat_file

    #Adding a progress bar (consider functionalizing to list only one time in future)
    # Calculate progress percentage
    progress_percent=$(( (curr_iter * 100) / max_iter ))

    # Build the progress bar
    bar_length=50
    completed_length=$(( progress_percent * bar_length / 100 ))
    remaining_length=$(( bar_length - completed_length ))
    completed_bar=$(printf "%${completed_length}s" | tr " " "#")
    remaining_bar=$(printf "%${remaining_length}s" | tr " " " ")

    # Print the progress bar
    printf "\rFasta file concatenation progress: [%s%s] %3d%%" "$completed_bar" "$remaining_bar" "$progress_percent"

    # Increment the counter
    curr_iter=$((curr_iter + 1))
done < $list_files
rm $list_files

cd ${here}

#Now we will run the command to generate the BLAST database
#For now we are using my conda environments, must change for 1.0 release
module purge
eval "$(conda shell.bash hook)"

conda activate NGStools

if [ "$fna_type" == true ]
then
    makeblastdb -in $tmp_dir$concat_file -dbtype nucl -out $fna_blastdb \
    -title nucl_db -max_file_sz 3GB
elif [ "$faa_type" == true ]
then
    makeblastdb -in $tmp_dir$concat_file -dbtype prot -out $faa_blastdb \
    -title prot_db -max_file_sz 3900000000
fi

#By default we will delete the contents of the tmp directory and just keep the BLAST
#db files in the reference data directory if needed for future use
if [ "$delete_tmp" == true ]
then
    printf "Removing temporary directory!\n\n"
    rm -rf $tmp_dir
fi
