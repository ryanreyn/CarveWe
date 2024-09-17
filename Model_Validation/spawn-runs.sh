#!/bin/bash

here=`pwd`
jobslist=${here}'/hq-matti-genomes.csv' # can change this to title of your joblist file
there='Runs'

#------Create run folder and move unique files into the run folder----------------------------------#
if [[ -d ${there} ]]
then
rm -r ${there}
mkdir ${there}

else
mkdir ${there}

fi

cp ${jobslist} ${there}

#-------------- Determine number of jobs you need to run  ----------------
#---------------------------------#
num_jobs=`wc -l ${jobslist} | awk '{print $1 }'`

echo 'Number of runs: ' ${num_jobs} '...'

#-------------- Loop over all polys (start after headed line)  ----------------------------#
iter=0

while [ ${iter} -lt ${num_jobs} ]
do
   let iter=${iter}+1
   curr_genome=`head -${iter} ${jobslist} | tail -1`
   runname="run-${curr_genome}"
   ####
   
   echo 'Run #'${iter}': creating files for '${runname}'...'
   
   #Create destination file names
   genPY=${there}"/run-${curr_genome}"'.py'
   genSLURM=${there}'/submit-'${curr_genome}'.slurm'
   
   #---- copy template and data files to runs folder 
cp -f ${here}'/batch-rescue-single-addition.py' ${genPY}
cp -f ${here}'/submit-batch-single-addition.slurm' ${genSLURM}
   # # # NOTE: This is currently making all files in the /Runs/ directory.
   # # # You can split the files to have them save in separate folders per job as well.
   
# ------- Editing the python code file ------- #
   joboutname=${runname}'_out'
   joberrname=${runname}'_err'
   
   sed -i -e s@sub_genome@${curr_genome}@g\
          -e s@sub_runname@${runname}@g    ${genPY}
  
            
# ------- Editing the Slurm files ------- #
joboutname=${there}'/'${runname}'_out'
joberrname=${there}'/'${runname}'_err'

sed -i -e s@sub_out@${joboutname}@g  \
       -e s@sub_err@${joberrname}@g  \
       -e s@sub_python@${genPY}@g    ${genSLURM}

echo 'Run files successfully modified'
  
done