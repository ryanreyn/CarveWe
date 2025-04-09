#!/bin/bash
#This is a preliminary functionality to be built out for CarveWe to enable a user to
#provide their own genome(s) and get the closest hit to a SOM cluster via BLAST
#Establish some flags for reading user input and providing help
fasta_file=""
is_dir="true"
faa_type="false"
fna_type="false"

while getopts "p:d:ano:t:b:" args; do
    case "${args}" in
        p) fasta_file=${OPTARG};;
        d) is_dir=${OPTARG};;
        a) faa_type='true';;
        n) fna_type='true';;
        o) out_dir=${OPTARG};;
        t) num_threads=${OPTARG};;
        b) blastdb_file=${OPTARG};;
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

outfile="${out_dir}/blast-to-carvewe.tsv"

#Run the appropriate blastn or blastp command depending on the n/a flag passed to
#the top level process call
if [ "$fna_type" == true ]
then
    echo "Running blastn"
    blastn -num_threads ${num_threads} \
    -query ${combined_fasta} \
    -db ${blastdb_file} \
    -outfmt 6 \
    -perc_identity 95 \
    -qcov_hsp_perc 100 \
    -num_alignments 1 \
    -out ${outfile} 2>/dev/null
elif [ "$faa_type" == true ]
then
    echo "Running blastp"
    blastp -num_threads ${num_threads} \
    -query ${combined_fasta} \
    -db ${blastdb_file} \
    -outfmt 6 \
    -perc_identity 95 \
    -qcov_hsp_perc 100 \
    -num_alignments 1 \
    -out ${outfile} 2>/dev/null
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
