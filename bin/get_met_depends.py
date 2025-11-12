#generate a dictionary with values in own media and media of 1st recipe, average media recipes at flux 1000,
#complete media recipes, own media recipes with flux set to 1000 and summed external flux after optimization in own media recipe
#look at peptides, inorganic, and amino acids

import os
import numpy as np
import pandas as pd
from scipy import stats
import optparse

from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.medium import minimal_medium

#Set up argument parsing
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-n', '--ensemble-size', type='int', dest='num_models',
                      help="provide the number of models to generate per ensemble")
    parser.add_option('-v', dest='verbose', action='store_true',
                      help="set this flag if you want to see all program output")
    parser.add_option('-w', '--work-dir', dest='work_dir',
                      help="provide the working directory that houses xml and reaction info files")
    parser.add_option('-o', '--media-dir', dest='media_dir',
                      help="provide the name of the directory you want the media output to be directed to")
    parser.add_option('-d', '--data-dir', dest='data_dir',
                      help="provide the name of the directory housing necessary data files")
    parser.add_option('-g', '--genome', dest='genome', default=None,
                      help="OPTIONAL: Process only files for this specific genome (for parallel mode). "
                           "If not provided, processes all genomes and creates consolidated output.")
    opts, args = parser.parse_args()
    #Convert input opts and args to python variables
    verbose=opts.verbose
    work_dir=opts.work_dir
    media_dir=opts.media_dir
    num_models=opts.num_models
    data_dir=opts.data_dir
    genome_filter=opts.genome


#step 1
#generate a list with all the genome names and one with all the ensemble sizes from the filenames in a folder

# PERFORMANCE FIX: Single-genome mode vs. batch mode
# If --genome is specified, only process that genome (parallel mode - no directory scanning overhead)
# If --genome is NOT specified, process all genomes (batch mode - original behavior)
if genome_filter:
    # Single-genome mode: No directory scanning, just use the provided genome name
    genomes = [genome_filter]
    if verbose:
        print(f"[SINGLE-GENOME MODE] Processing only genome: {genome_filter}")
else:
    # Batch mode: Scan directory for all genomes (original behavior)
    #list of directory with all the rxn info files
    #rxn info files have ensemble states and rxn names
    directory = os.fsencode('%s/xml_files/' %(work_dir))
    os.listdir(directory)

    #read in an individual file
    #loop through all the files in the folder
    genomes = []
    for file in os.listdir(directory):
      filename = os.fsdecode(file)
      if filename.endswith('.xml'):

        #parse name of model
        name = filename.replace('.xml', '') #filename without .xml

        if name not in genomes:
          genomes.append(name)

    if verbose:
        print(f"[BATCH MODE] Processing {len(genomes)} genomes from directory scan")



def restore_rxn_states(model, tup_list):
  #restores the reaction states from a list of tuples
  i = 0 #indexing variable
  for x in model.reactions:
    tup = tup_list[i]
    x.lower_bound = tup[0]
    x.upper_bound = tup[1]
    i += 1

#a function that turns off all nutrient imputs for a given model
def no_nutrients(model):
    for x in model.reactions:
        if x.id.startswith('EX') == True:
            model.reactions.get_by_id(x.id).lower_bound = 0.

#function to set media based on the flux values of a dataframe column:
#dataframe must have indexes with the external rxns on them
def column_to_media(model, df, column):
  no_nutrients(model) # sets all external reatction fluxes to 0
  for i in range(0, len(df.index)): #sets external fluxes based on the media df
    if np.isnan(df.iloc[i, column]) == False:
      try:
          model.reactions.get_by_id(df.index[i]).lower_bound = -1 * df.iloc[i, column]
      except KeyError:
          continue

#a function that turns on all nutrient imputs for a given model
def all_nutrients(model):
    for x in model.reactions:
        if x.id.startswith('EX') == True:
            model.reactions.get_by_id(x.id).lower_bound = -1000.

#function to set media based to 1000 based on the presence of a metabolite in a dataframe column:
#dataframe must have indexes with the external rxns on them
def column_to_media_1000(model, df, column):
  no_nutrients(model) # sets all external reatction fluxes to 0
  for i in range(0, len(df.index)): #sets external fluxes based on the media df
    if np.isnan(df.iloc[i, column]) == False:
      try:
          model.reactions.get_by_id(df.index[i]).lower_bound = -1000
      except KeyError:
          continue

#function to set media based to 1000 based on the presence of a metabolite in a dataframe column:
#dataframe must have indexes with the external rxns on them
def column_to_media_half_flux(model, df, column):
  no_nutrients(model) # sets all external reatction fluxes to 0
  for i in range(0, len(df.index)): #sets external fluxes based on the media df
    if np.isnan(df.iloc[i, column]) == False:
      try:
          model.reactions.get_by_id(df.index[i]).lower_bound = -1 * df.iloc[i, column]/2
      except KeyError:
          continue

#function to set media based on the flux values of a dataframe column:
#flux values will be set to 1/2 original values for metabolites in given category
#otherwise they will be set to -1000
#dataframe must have indexes with the external rxns on them
def category_media(model, df, column, cat_list):
  no_nutrients(model) # sets all external reatction fluxes to 0
  for i in range(0, len(df.index)): #sets external fluxes based on the media df
    if np.isnan(df.iloc[i, column]) == False:
      #if metabolite is in a category, set the flux to half the original input flux
      if df.index[i] in cat_list:
        try:
            model.reactions.get_by_id(df.index[i]).lower_bound = -1 * df.iloc[i, column]/2
        except KeyError:
            continue
      #otherwise set the value to -1000
      else:
        try:
            model.reactions.get_by_id(df.index[i]).lower_bound = -1000
        except KeyError:
            continue

growth_info = {} #initialize dictionary
sensitivity_info = {} #initialize second dictionary
cats = ['Carboxylic Acid', 'Other', 'Ketones/Aldehydes',
       'Amino Acids/Derivatives',
       'Nucleobases/Nucleosides/Nucleotides/Derivatives', 'Peptides',
       'Inorganic', 'Organic Sulfur', 'Alcohol', 'Amines/Amides',
       'B Vitamins', 'Carbohydrates/Derivatives',
       'Phospholipids/Fatty Acids/Triglycerides'] #sets the metabolite categories

# PERFORMANCE FIX: Read genome-specific input files in single-genome mode
# Determine input file suffix based on mode (must match convert-media-output.py output)
if genome_filter:
    # Single-genome mode: Read genome-specific files created by convert-media-output.py
    file_suffix = f"_{genome_filter}"
    if verbose:
        print(f"[SINGLE-GENOME MODE] Reading genome-specific input files with suffix: {file_suffix}")
else:
    # Batch mode: Read consolidated files (original behavior)
    file_suffix = ""
    if verbose:
        print(f"[BATCH MODE] Reading consolidated input files")

#pull in the media data with headers
data = pd.read_csv('%s/reformatted_media_data%s.csv' %(media_dir, file_suffix), index_col = 0)
#data = data.fillna(0)

#pull in the classified metabolite data
met_class = pd.read_csv('%s/classified_metabolites.csv' %(data_dir), index_col=0)

#load in the metabolite names dictionary
#pull in the data
met_names_df = pd.read_csv('%s/all_max_flux_met_names%s.csv' %(media_dir, file_suffix), index_col = 0)
met_names = met_names_df.squeeze()

#loop though the genomes
for genome in genomes:
#print(genome)
  #look at the genomes in that SOM category
  if genome in data['model_name'].to_list():

    growth_info[genome] = {}
    growth_info[genome]['growth_rates_own'] = [] #growth rate in own media

    count = 0

    #also set up a dataframe to store the sensitivity values calculated using Eq 1
    sensitivity_info[genome] = {}

    #import model and reaction states
    model = read_sbml_model('%s/xml_files/%s.xml' %(work_dir, genome))
    rxn_states = pd.read_csv('%s/ensemble_rxn-info_files/%s_rxn-info.csv' %(work_dir, genome), header = None, index_col=0)

    #subset temp dataframe for only the single genome
    min = data[data['model_name'] == genome].drop(columns= ['model_name', 'model_number', 'media_type', 'ensemble_size']).T

    # print(min)

    #store original bounds for fluxes for this model:
    original_bounds = []
    for x in model.reactions:
      lb = x.lower_bound
      ub = x.upper_bound
      original_bounds.append((lb, ub))

    #loop through all the individual models
    while count < num_models:

      #restore original bounds
      restore_rxn_states(model, original_bounds)

      #set reaction bounds to 0 for reactions that have reaction state of 0
      i = 0
      for x in model.reactions:
        is_present = rxn_states.iloc[i,count]
        if is_present == 0:
          x.lower_bound = 0.
          x.upper_bound = 0.
        i +=1


      #own media
      #set the fluxes of the media components based on the recipe:
      column_to_media(model, min, count)
      #slim optimize and record value:
      solution = model.slim_optimize()
      growth_info[genome]['growth_rates_own'].append(solution)


      for met_cat in cats:
        if met_cat not in growth_info[genome].keys(): #check if the key exists for this category and if not initialize
          growth_info[genome][met_cat] = [] #growth rates for a specific category, with 1/2 fluxes for that category and -1000 for everything else
          sensitivity_info[genome][met_cat] = []
        #media based on a category
        cat_list = met_names[met_names.isin(met_class[met_class["Higher Level Classification"] == met_cat].index)].index.to_list()
        #set fluxes based on a column and the given category:
        category_media(model, min, count, cat_list)
        #slim optimize and record value:
        solution = model.slim_optimize()
        growth_info[genome][met_cat].append(solution)

        #record the sensitivity value of each met cat as well
        sens_value = 2*(1-solution/growth_info[genome]['growth_rates_own'][count])
        if sens_value > 1:
          sens_value = 1
          sensitivity_info[genome][met_cat].append(sens_value)
        else:
          sensitivity_info[genome][met_cat].append(sens_value)


      count += 1

print(growth_info)

# PERFORMANCE FIX: Genome-specific output files in single-genome mode to prevent overwrites
# Output file suffix matches the mode used throughout
if genome_filter:
    # Single-genome mode: Add genome name to output files (prevents parallel job overwrites)
    output_suffix = f"_{genome_filter}"
    if verbose:
        print(f"[SINGLE-GENOME MODE] Writing genome-specific output files with suffix: {output_suffix}")
else:
    # Batch mode: Use standard filenames (original behavior)
    output_suffix = ""
    if verbose:
        print(f"[BATCH MODE] Writing consolidated output files")

cat_growth_info_df = pd.DataFrame.from_dict(growth_info)
cat_growth_info_df.to_csv('%s/model_growth_rate_info_by_met_category_all_metabolites%s.csv' %(media_dir, output_suffix))

#output sensitivity info to a .csv as well
sens_growth_info_df = pd.DataFrame.from_dict(sensitivity_info)
sens_growth_info_df.to_csv('%s/model_sensitivity_by_met_category%s.csv' %(media_dir, output_suffix))
