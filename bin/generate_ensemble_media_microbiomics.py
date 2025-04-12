# this version is edited to run on the cluster, it's then coppied and pasted into a text file
# this version is edited to run generally for a single file that has a genome name input as an agument
# run like this: 'python thisscript.py genome_name'

import numpy as np
import pandas as pd
import optparse

from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.medium import minimal_medium
import os 
import sys
from datetime import datetime

#Set up argument parsing
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-n', '--ensemble-size', type='int', dest='num_models',
                      help="provide the number of models to generate per ensemble")
    parser.add_option('-v', dest='verbose', action='store_true',
                      help="set this flag if you want to see all program output")
    parser.add_option('-d', '--work-dir', dest='work_dir',
                      help="provide the working directory that houses xml and reaction info files")
    parser.add_option('-o', '--media-dir', dest='media_dir',
                      help="provide the name of the directory you want the media output to be directed to")
    opts, args = parser.parse_args()
    #Convert input opts and args to python variables
    num_models=opts.num_models
    verbose=opts.verbose
    work_dir=opts.work_dir
    media_dir=opts.media_dir
    genome = args[0]


#define a function that outputs a dataframe with all minimal media options
def find_min_medias(model):
  #optimize the model
  solution = model.slim_optimize()

  #run the minimal media function to create an ensemble of up to 10 minimal media outputs
  min = minimal_medium(model, solution, minimize_components=10)

  #this minimal media function maximizes flux but does not minimize number of components
  min1 = minimal_medium(model, solution, minimize_components=False)

  #combine min and min1 into a dataframe
  min1_df = min1.to_frame()
  min1_df.rename(columns={0:'max_flux'}, inplace=True)
  min_df = pd.concat([min1_df, min], axis=1)
  min_df

  return min_df

#this functin adds metabolite names to a dataframe with metabolite ids as the index
def add_met_names(model, min_df):
  #add metabolite names in human readable format
  min_df_names = min_df.copy(deep=True)
  for index in min_df_names.index:
    met = index.replace("EX_", "")
    name = model.metabolites.get_by_id(met).name
    min_df_names.loc[index, "Met_Name"] = [name]

  return min_df_names

def restore_rxn_states(model, tup_list):
  #restores the reaction states from a list of tuples
  i = 0 #indexing variable
  for x in model.reactions:
    tup = tup_list[i]
    x.lower_bound = tup[0]
    x.upper_bound = tup[1]
    i += 1 


#We need to import the .xml file or SBML model as well as the reaction states file
model = read_sbml_model('%sxml_files/%s.xml' %(work_dir, genome))
rxn_states = pd.read_csv('%sensemble_rxn-info_files/%s_rxn-info.csv' %(work_dir, genome), header = None, index_col=0)

headers = pd.DataFrame()
min_media = pd.DataFrame()

original_bounds = []
for x in model.reactions:
  lb = x.lower_bound
  ub = x.upper_bound
  original_bounds.append((lb, ub))

for j in range(0, num_models): #column indexing variable
  i = 0 #reaction indexing variable

  #restore original bounds
  restore_rxn_states(model, original_bounds)

  #set reaction bounds to 0 for reactions that have reaction state of 0
  for x in model.reactions:
    is_present = rxn_states.iloc[i,j]
    if is_present == 0:
      x.lower_bound = 0.
      x.upper_bound = 0.
    i +=1

  #generate minimal media recipe
  try:
    min = find_min_medias(model)
    min_media = pd.concat([min_media, min], axis = 1)

    #generate a header to save the metadata
    model_name = []
    model_number = []
    media_type = []
    ensemble_size = []
    name = genome
    for column in min.columns:
      c = 0
      if column == "max_flux":
        c = 1
      media_type.append(c)
      model_name.append(name)
      model_number.append(j)
      ensemble_size.append(num_models)
    header = pd.DataFrame([model_name, model_number, media_type, ensemble_size], columns = min.columns, index = ['model_name', 'model_number', 'media_type', 'ensemble_size'])
  
  except AttributeError:
    print('Media could not be found for %s ensemble' %(genome))

  #append header info to headers dataframe
  headers = pd.concat([headers, header], axis=1)

#add metabolite names
all_media = add_met_names(model, min_media)

all_media.to_csv("%s%s/%s_ensemble_media.csv" %(work_dir, media_dir, genome))
headers.to_csv("%s%s/%s_ensemble_header.csv" %(work_dir, media_dir, genome))
