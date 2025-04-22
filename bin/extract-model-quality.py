import numpy as np
import pandas as pd
import optparse

#Set up argument parsing
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file-list', type='int', dest='file_list',
                      help="provide the file that lists all of the genome names to run the program on")
	#Convert input opts and args to python variables
    file=opts.file_list

#This script will take files from a file list and extract their model consensus
#Then we'll add it to a csv file for analysis
filelist=open(file,"r")

save_data=[["Genome","Consensus"]]
for curr_file in filelist:
	curr_file=curr_file.strip()
	genome=curr_file.replace("_rxn-info.csv","")
	curr_data=pd.read_csv(curr_file)
	per_rxn_consensus=curr_data.sum(numeric_only=True,axis=1)/60
	avg_consensus=per_rxn_consensus.mean()
	data_row=[genome,avg_consensus]
	save_data.append(data_row)

save_df=pd.DataFrame(save_data)
pd.DataFrame.to_csv(save_df,"quality-data.csv",header=False,index=False,quoting=None)
