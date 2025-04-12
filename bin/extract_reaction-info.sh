#!/bin/bash
#This code is designed to take a .xml ensemble file and extract the reaction IDs 
#and associated presence/absence states for all models in the ensemble
printf "Running xml extraction:\n\n"

xml_dir="xml_files"

while getopts "o:x:" args; do
    case "${args}" in
        o) out_dir=${OPTARG};;
        x) xml_dir=${OPTARG};;
    esac
done

#Establish the name of the files list and create that file and a directory to write 
#all of the output files into
filelist=${out_dir}"filenames.txt"
dirname=${out_dir}"ensemble_rxn-info_files"

if [[ -d $dirname ]]
then
    rm -r $dirname
    mkdir $dirname

elif [[ ! -d $dirname ]]
then
    mkdir $dirname
    touch $filelist

fi

#We will now construct the files list from all .xml files in the current directory
ls ${xml_dir} > $filelist

#We are going to do this in an iterative approach so it can be done on all of the ensembles
while read file
do
    #Establish a save file name for the reaction info and the names of the temporary files
    #for each piece to save
    printf "Extracting reaction info from $file into $dirname\n"

    prefix=`echo $file | cut -d . -f 1 | cut -d \/ -f 2`
    savename=$dirname"/"$prefix"_rxn-info.csv"
    reactions="reaction_IDs_tmp.tsv"
    states="ensemble_states_tmp.csv"

    #First search for the reaction IDs and save it in a temporary file
    `grep -e "reaction metaid" $xml_dir/$file | cut -d \" -f 2 > $reactions`

    #Now we will extract the ensemble states associated with each reaction ID and save as 
    #a temporary file
    `grep -e "ENSEMBLE_STATE" $xml_dir/$file | cut -d : -f 2 | cut -d '<' -f 1 | sed 's/ /,/2g' > $states`

    #Using the paste command to combine the two files and then removing our temporary files
    `paste -d , $reactions $states > $savename`

    #Remove spaces from the files
    `sed -i 's/\s//g' $savename`

    rm $reactions
    rm $states

done < $filelist

#Remove any temporary files
rm $filelist
