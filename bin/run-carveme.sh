#!/bin/bash
#This script is designed to take user input from the higher level CarveWe bash script
#and perform several steps of the CarveWe metabolic model annotation pipeline
fna_type="false"
faa_type="false"
max_processes=1
ensemble_size=2

printf "Running CarveMe\n\n"

while getopts "i:o:e:anp:" args; do
    case "${args}" in
        i) fasta_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        e) ensemble_size=${OPTARG};;
        n) fna_type='true';;
        a) faa_type='true';;
        p) max_processes=$((OPTARG - 1));;
    esac
done

#Activating the CarveMe conda environment **will need to change this for final release**
module purge

eval "$(conda shell.bash hook)"

conda activate /project/nmherrer_110/tools/.conda/envs/CarveMe

#Setting up the destination directory for the xml files
xml_dir=$out_dir"/xml_files"

if [ -d $xml_dir ]
then
    rm -r $xml_dir
    mkdir $xml_dir
else
    mkdir $xml_dir
fi


#List out all of the fasta files to a temporary file to draw from to run CarveMe
tmp_filelist=${out_dir}"/tmp-fasta-list.txt"
ls ${fasta_dir}*.f[an]a > $tmp_filelist

#Set up some counting variables for multiprocessing and progress display
current_processes=0
curr_iter=1
max_iter=`wc -l $tmp_filelist | cut -d ' ' -f 1`

# Initialize the progress bar
progress_percent=$(( (0 * 100) / max_iter ))
bar_length=50
completed_length=$(( progress_percent * bar_length / 100 ))
remaining_length=$(( bar_length - completed_length ))
completed_bar=$(printf "%${completed_length}s" | tr " " "#")
remaining_bar=$(printf "%${remaining_length}s" | tr " " " ")
printf "\rModel building in progress: [%s%s] %3d%%" "$completed_bar" "$remaining_bar" "$progress_percent"
printf "\n"

#Actual command loop to define fasta and xml files and execute CarveMe
while read fasta_file
do
    while [[ $current_processes -ge $max_processes ]]
    do
        wait -n # Waiting for any process to complete to add more processes
        current_processes=$((current_processes - 1))
    done

    xml_file=`echo $fasta_file | sed "s/\.f[an]a/\.xml/g; s/.*\///g"`

    #printf "Creating CarveMe ensemble ${xml_file} in ${xml_dir} with ${ensemble_size} models\n\n"
    
    if [ "$fna_type" == 'true' ]
    then
        (
        carve -v --dna -o "${xml_dir}/${xml_file}" -n $ensemble_size $fasta_file > /dev/null

        #Adding a progress bar instead of letting standard output hit terminal
        # Calculate progress percentage
        progress_percent=$(( (curr_iter * 100) / max_iter ))

        # Build the progress bar
        bar_length=50
        completed_length=$(( progress_percent * bar_length / 100 ))
        remaining_length=$(( bar_length - completed_length ))
        completed_bar=$(printf "%${completed_length}s" | tr " " "#")
        remaining_bar=$(printf "%${remaining_length}s" | tr " " " ")

        # Print the progress bar
        printf "\rModel building in progress: [%s%s] %3d%%" "$completed_bar" "$remaining_bar" "$progress_percent"
        printf "\n"
        # Increment the counter
        curr_iter=$((curr_iter + 1))
        ) &
        ((current_processes ++ ))
    else
        (
        carve -v -o "${xml_dir}/${xml_file}" -n $ensemble_size $fasta_file \
        > /dev/null
        #printf "\n"

        #Adding a progress bar instead of letting standard output hit terminal
        # Calculate progress percentage
        completed_files=`ls $xml_dir | wc -l | cut -d ' ' -f 1`
        progress_percent=$(( (completed_files * 100) / max_iter ))

        # Build the progress bar
        bar_length=50
        completed_length=$(( progress_percent * bar_length / 100 ))
        remaining_length=$(( bar_length - completed_length ))
        completed_bar=$(printf "%${completed_length}s" | tr " " "#")
        remaining_bar=$(printf "%${remaining_length}s" | tr " " " ")

        # Print the progress bar
        printf "\rModel building in progress: [%s%s] %3d%%" "$completed_bar" "$remaining_bar" "$progress_percent"
        printf "\n"
        ) &
        ((current_processes ++ ))
    fi
done < $tmp_filelist

#Make sure to wait for any processes that are still running to finish
while [[ $current_processes -gt 0 ]]
do
    wait -n
    current_processes=$((current_processes - 1))
done

printf "\n"

#Now we will process the output .xml files with our custom JSON parsing file to
#create our ensemble reaction info files
#extract_reaction-info.sh -x $xml_dir -o $out_dir
python "${scripts_dir}/extract_reaction-info.py" -x $xml_dir -o $out_dir

#Remove any temporary files made in the course of the subprocess
rm $tmp_filelist

#Deactivate conda environment before releasing back to main process
conda deactivate
