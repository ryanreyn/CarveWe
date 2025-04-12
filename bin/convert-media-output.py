#This script will take the raw media output and isolate the max flux recipes for each
#ensemble model for each genome, then reformat the data
import numpy as np
import pandas as pd
import os
import sys
import optparse

#Set up an options parser to pull in necessary info from higher level process
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-v', dest='verbose', action='store_true',
                      help="set this flag if you want to see all program output")
    parser.add_option('-d', '--work-dir', dest='work_dir',
                      help="provide the working directory that houses xml and reaction info files")
    parser.add_option('-o', '--media-dir', dest='media_dir',
                      help="provide the name of the directory you want the media output to be directed to")
    opts, args = parser.parse_args()
    #Convert input opts and args to python variables
    verbose=opts.verbose
    work_dir=opts.work_dir
    media_dir=opts.media_dir

#step 1
#generate a list with all the genome names and one with all the ensemble sizes

#list of directory with all the rxn info files
#rxn info files have ensemble states and rxn names
directory = os.fsencode('%s' %(media_dir))
os.listdir(directory)

#read in an individual file
#loop through all the files in the folder
genomes = []
for file in os.listdir(directory):
  filename = os.fsdecode(file)

  #parse name of model
  name = filename.replace('.csv', '') #filename without .csv
  split = name.split('_')
  name = '_'.join(split[:-1])
  split = split[:-2] #gets genome name without 'ensemble'
  strname = '_'.join(split) #puts ensemble name into underscore separated format

  if strname not in genomes:
    genomes.append(strname)


len(genomes)


#step 2
#generate unsubseted dataframes of maxflux metabolites without averaging

#Generate dictionaries for rank abundance curves:
#dictionaries will contain: name, count, count per genome
#total number of media recipes per genome and total media recipes will be stored separately
#it will also generate a dataframe containing the average fluxes of the core reactions for each enseble, with metabolites as rows and genomes as columns

met_info = {} #nested dictionary with information about each metabolite
genome_totals = {} #dictionary that contains all the total number of media recipes for each genome
total_recipes = 0 #variable that stores total number of media recipes
fluxes_df = pd.DataFrame() #dataframe with all the average fluxes for core reactions of each model
core = 0 #sets the cutoff for what is considered a core reaction, represented as a decimal between 0 and 1
max_flux = True #only use max_flux reactions
headers = pd.DataFrame()
met_names = {}

#initialize genome totals dictionary
for genome in genomes:
  genome_totals[genome]=0

for genome in genomes:
    #load files in order of ensemble size
    rxn_states = pd.read_csv('%s/%s_ensemble_media.csv' %(media_dir, genome), index_col=0)
    header = pd.read_csv('%s/%s_ensemble_header.csv' %(media_dir, genome), index_col=0)
    #rxn_states = rxn_states.fillna(0) #fill NaNs with 0's

    # if metabolite is not in dictionary, intialize the entry for that rxn state
    for index in rxn_states.index:
      if index not in met_info.keys():
        met_info[index] = {'name': rxn_states.loc[index, 'Met_Name'], 'count': 0, 'genome_count':{}}
        for genome in genomes:
          met_info[index]['genome_count'][genome] = 0

    # add metabolite names to a dictionary
    for i in rxn_states.index:
      if i not in met_names.keys():
        met_names[i] = rxn_states.loc[i, 'Met_Name']

    #if looking for max_flux only filter for those columns
    if max_flux == True:
      rxn_states = rxn_states.filter(regex='max_flux')
      header = header.filter(regex='max_flux')

    shape = rxn_states.shape #save size of dataframe, this will be used to count the number of media recipes

    #add the number of recipes to both the overal total and genome total
    genome_totals[genome] += shape[1]
    total_recipes += shape[1]

    #count the number of recipes containing each metabolite into a new column
    #rxn_states['total'] = rxn_states.loc[:, ~rxn_states.columns.isin(['Met_Name', 'mean'])].select_dtypes(np.number).gt(0).sum(axis=1)

    #calculate the percentage of recipes with each metabolite
    #rxn_states['percentage'] = rxn_states['total']/shape[1]

    #creates a temporary dataframe with mean flux values for core reactions
    #temp = rxn_states.query('percentage>%f' %core)['mean']
    #temp.rename(j, inplace = True)

    #adds the core rxn means to the avg_flux dataframe
    fluxes_df = pd.concat([fluxes_df, rxn_states], axis=1)

    #add headers to a file
    headers = pd.concat([headers, header], axis=1)

    #add metabolite counts for total and individual genomes
    #for index in rxn_states.index:
    #  met_info[index]['count'] += rxn_states.loc[index, 'total']
    #  met_info[index]['genome_count'][j] += rxn_states.loc[index, 'total']

#remove NaNs from the fluxes dataframe:
fluxes_df_filled = fluxes_df.fillna(0)
fluxes_df_filled.info()

#store all max flux dataframes
fluxes_df_filled.to_csv('%s/all_max_flux_diamond.csv' %(work_dir))
headers.to_csv('%s/all_max_flux_headers_diamond.csv' %(work_dir))

names = pd.DataFrame.from_dict(met_names, orient='index')
names.to_csv('%s/all_max_flux_met_names.csv' %(work_dir))

# pull in values and headers
max_flux = pd.read_csv('%s/all_max_flux_diamond.csv' %(work_dir), index_col = 0)
headers = pd.read_csv('%s/all_max_flux_headers_diamond.csv' %(work_dir), index_col = 0)

#replace 0s with NaNs in the metabolite dataframe
max_flux.replace(0, np.nan, inplace=True)

#replace column names with the model names
columns = headers.columns
headers.set_axis(headers.loc['model_name'].to_list(), axis='columns', inplace=True) 
#set column names to model names

#set back the index to the original index
headers.set_axis(columns, axis ='columns', inplace = True)

#concatinate the headers and the data
data = pd.concat([headers, max_flux])
data = data.T

data.to_csv('%s/reformatted_media_data.csv' %(work_dir))
