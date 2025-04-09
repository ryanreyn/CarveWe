#!/bin/bash
#This script is designed to take user input from the higher level CarveWe bash script
#and perform several steps of the CarveWe metabolic model annotation pipeline
printf "Running CarveMe\n\n"

while getopts "i:o:e:" args; do
    case "${args}" in
        i) fasta_dir=${OPTARG};;
        o) out_dir=${OPTARG};;
        e) ensemble_size=${OPTARG};;
    esac
done

#Activating the CarveMe conda environment **will need to change this for final release**
module purge

eval "$(conda shell.bash hook)"

conda activate /project/nmherrer_110/tools/.conda/envs/CarveMe

#We will extract the genome size (number of lines in eggnog file) and BiGG genes and pipe them to the output file regardless of if the carveme command runs or fails
#genome=`echo sub_egg_file | cut -d \- -f 1`
#Size=`sed '1d' sub_egg_file | wc -l`
#BiGG=`cut -f 17 sub_egg_file | sed '/^\s*$/d' | sed '1d' | wc -l`
#echo "Genome Size (number of genes): " $Size >> sub_output
#echo "Number of BiGG genes: " $BiGG >> sub_output
xml_dir=$out_dir"xml_files"

#List out all of the fasta files to a temporary file to draw from to run CarveMe
tmp_filelist=${out_dir}"tmp-fasta-list.txt"
ls ${fasta_dir}*.fna > $tmp_filelist

#Actual command loop to define fasta and xml files and execute CarveMe
while read fasta_file
do
    xml_file=`echo $fasta_file | sed "s/\.f[an]a/\.xml/g; s/.*\///g"`
    printf "Creating CarveMe ensemble ${xml_file} in ${xml_dir} with ${ensemble_size} models\n\n"
    carve -v --dna -o "${xml_dir}/${xml_file}" -n $ensemble_size $fasta_file
    printf "\n"
done < $tmp_filelist

#Now we will process the output .xml files with our custom JSON parsing file to
#create our ensemble reaction info files
extract_reaction-info.sh -x $xml_dir -o $out_dir

#Remove any temporary files made in the course of the subprocess
rm $tmp_filelist

#Deactivate conda environment before releasing back to main process
conda deactivate
