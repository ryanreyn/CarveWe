#!/bin/bash
#This is a small bash script for compiling batch output files
#Create empty file with desired name for compiled output
here=`pwd`
there=${here}"/Runs"
jobslist=${here}"/hq-matti-genomes.csv"

recovery_file="../Output/batch-predicted_media_single_add_test_recovery.csv"
rescue_file="../Output/batch-predicted_media_single_add_test_rescued_metabolies.csv"

if [[ -f ${recovery_file} ]] || [[ -f ${rescue_file} ]]
then
    rm ${recovery_file}
    rm ${rescue_file}

    touch ${recovery_file}
    touch ${rescue_file}

else
    touch ${recovery_file}
    touch ${rescue_file}

fi

while read line
do
    curr_recovery=${there}"/run-${line}_recovery.csv"
    curr_rescue=${there}"/run-${line}_rescued_metabolites.csv"

    cat ${curr_recovery} >> ${recovery_file}
    cat ${curr_rescue} >> ${rescue_file}

done < ${jobslist}