import numpy as np
import pandas as pd
import optparse
import os
import seaborn as sns

#Set up argument parsing
if __name__ == '__main__':
    parser = optparse.OptionParser()
	parser.add_option('-d', '--directory', type='int', dest='file_dir',
					  help="provide the directory that contains the model files")
	parser.add_option('-w', '--work-directory', type='int', dest='work_dir',
					  help="provide the working directory where the compiled output should go")
	#Convert input opts and args to python variables
    file=opts.file_list
	file_dir=opts.file_dir
	work_dir=opts.out_dir

#This script will take files from a directory and extract their model consensus
#Then we'll add it to a csv file for analysis
directory = os.fsencode('%s/%s' %(work_dir, file_dir))
os.listdir(directory)


#read in an individual file
#loop through all the files in the folder
for file in os.listdir(directory):
  filename = os.fsdecode(file)

save_data=[["Genome","Consensus"]]
for file in os.listdir(directory):
	curr_file=os.fsdecode(file)
	if filename.endswith('_rxn-info.csv'):
    	#parse name of model
    	genome=filename.replace('_rxn-info.csv', '') #filename without ending

	genome=curr_file.replace("_rxn-info.csv","")
	curr_data=pd.read_csv(curr_file)
	per_rxn_consensus=curr_data.sum(numeric_only=True,axis=1)/60
	avg_consensus=per_rxn_consensus.mean()
	data_row=[genome,avg_consensus]
	save_data.append(data_row)

save_df=pd.DataFrame(save_data)
pd.DataFrame.to_csv(save_df,"%s/quality-data.csv" %(work_dir),header=False,index=False,quoting=None)

#We want to extend this script to generate some seaborn plots of metrics related to
#model quality (potential for adding an exploratory tool for picking ensemble size)
quality_df=save_df

quality_hist=sns.histplot(data=quality_df, x="Consensus", y="Count", kde=True)
quality_hist.savefig("model-quality-histogram.png")