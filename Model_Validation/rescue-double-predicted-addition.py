#this version finds all the unique compounds in the predicted media recipes, adds them one at a time
#and then solves for all of the matti compounds
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.medium import minimal_medium


import numpy as np
import pandas as pd


#import matplotlib.patches as pat
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from matplotlib import colors


#from matplotlib import cm
#from PIL import Image, ImageDraw
import os
import fnmatch
from pandas.io.parsers.python_parser import count_empty_vals

import ast

import re

#Try to suppress the cobra.medium log messages that are bloating the output
import logging
cobra_logger = logging.getLogger('cobra.medium')
cobra_logger.setLevel(logging.CRITICAL)

#Output a print statement
print("Successfully loaded packages!")

#name of this run/run type:
run_name = "predicted_media_double_add_test"

#step 1
#generate a list with all the genome names and one with all the ensemble sizes

#list of directory with all the rxn info files
#rxn info files have ensemble states and rxn names
directory = os.fsencode('/project/nmherrer_110/acweiss/ensemble_media/matti_media')
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

#store reaction states in a dictionary
rxn_info = {}
for genome in genomes:
  rxn_states = pd.read_csv('/project/nmherrer_110/acweiss/ensemble_media/matti-data_ensemble_rxn-info_files/%s_rxn-info.csv' %genome, header = None, index_col=0)
  mean_freq = rxn_states.mean(axis='columns').mean(axis='index') #calculates the mean reaction frequency
  stdev = rxn_states.mean(axis='columns').std(axis='index') #standard deviation of reaction frequencies
  num_rxns = rxn_states.shape[0] #the total number of reactions for the genome
  median_freq = rxn_states.mean(axis='columns').median(axis='index') #calculates the median reaction frequency
  #calcuate the number of core reactions
  core = .8 #threshold for core reactions
  freq = rxn_states.mean(axis='columns') # save the mean reaction frequencies in a series
  core_count = freq[freq >= core].count() # count the number above the core threshold value
  rxn_info[genome] = {'mean_freq': mean_freq, 'stdev' : stdev, 'num_rxns': num_rxns, 'median_freq': median_freq, 'core_rxns': core_count}

#create a dataframe from the dictionary
rxn_df = pd.DataFrame.from_dict(rxn_info, orient = "index")

len(rxn_df[(rxn_df.mean_freq >= 0.8)])

#load in all the dependancy files
#load matti vitamins
vitamins = pd.read_csv('../Data/matti_vitamins.csv', index_col = 0)
#load the reaction information dataframe
rxn_freqs = pd.read_csv('../Data/matti_rxn_info.csv', index_col = 0)
#load in the lab IDs file
lab_mets = pd.read_csv('../Data/matti_lab_metabolite_ids.csv', index_col = 0)
#import failing genomes file
failing_genomes = pd.read_csv('../Data/matti_rescues_failing_genomes_20241010.csv', index_col = 0)


#set high quality genomes based on index of failing genomes file
hq_genomes = failing_genomes.index

#Set all the carbon containing reactions to 0 and turn on all the other reactions


#set path locations:
xml_path = '/project/nmherrer_110/acweiss/ensemble_media/matti-data_xml-files'
data_path = '/project/nmherrer_110/acweiss/ensemble_media/matti_media'
rxn_path = '/project/nmherrer_110/acweiss/ensemble_media/matti-data_ensemble_rxn-info_files'


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

# Function to identify carbon-containing metabolites
def contains_carbon(metabolite):
    if metabolite.formula:
        # Regex to match 'C' not followed by a lowercase letter
        if re.search(r'C(?![a-z])', metabolite.formula):
            return True
    return False

#function that finds the unique compounds in the predicted media
#ie the carbon compounds not accounted for in vitamins or mattis compounds
def find_media_rxns(genome, model, data_path):
  #import model and reaction states
  media = pd.read_csv('%s/%s_ensemble_media.csv'%(data_path, genome), index_col = 0)

  #first we need to filter for the metabolites that are are in the max_flux recipes only
  media = media.filter(like="max_flux").dropna(how="all")
  #then we filter for metabolite presence in at least 50% of recipes
  media = media.dropna(thresh=int(media.shape[1] * 0.5))

  #next we need to filter out anything that's not a carbon source
  #generate a list of all the non C sources
  noC = []
  for metid in media.index:
    for metabolite in model.reactions.get_by_id(metid).metabolites:
        if contains_carbon(metabolite) == False:
            noC.append(metid)
  #remove those rows
  media = media.drop(index=noC)

  #now try dropping the vitamins
  media = media.drop(index=vitamins.index, errors='ignore')

  #in this case don't drop the matti compounds
  #now try droppint the matti compounds:
  #media = media.drop(index=lab_mets['Exchange_rxn'].values, errors='ignore')

  return media

#initialize dictionary
growth_info = {} #initialize dictionary

for genome in hq_genomes:
    #import model and reaction states
  model = read_sbml_model('%s/%s.xml'%(xml_path, genome))
  rxn_states = pd.read_csv('%s/%s_rxn-info.csv' %(rxn_path, genome), header = None, index_col=0)

  #initialize for genome
  growth_info[genome] = {}
  #growth_info[genome]['growth_rates_no_C'] = [] #growth rate in own media

  #store original bounds for fluxes for this model:
  original_bounds = []
  for x in model.reactions:
    lb = x.lower_bound
    ub = x.upper_bound
    original_bounds.append((lb, ub))

  #generate dataframe of predicted media reactions which are unique to that recipe
  #meaning they are carbon containing but not matti compounds or vitamins
  pred_med = find_media_rxns(genome, model, data_path)

  #set rescue compound based on best rescue
  rescue_comp = failing_genomes.loc[genome, "best_rescue"]

  #loop through for each predicted media component
  for predrxn in pred_med.index:
    growth_info[genome][predrxn] = {}

    count = 0

    #only use rxns other than the rescue comp:
    if predrxn != rescue_comp:

      #loop through all the individual models
      while count < 60:
                #restore original bounds
        restore_rxn_states(model, original_bounds)

        #initialize for model
        growth_info[genome][predrxn][count] = {}

        #set reaction bounds to 0 for reactions that have reaction state of 0
        i = 0
        for x in model.reactions:
          is_present = rxn_states.iloc[i,count]
          if is_present == 0:
            x.lower_bound = 0.
            x.upper_bound = 0.
          i +=1

        #save positive control, all rxns on
        solution = model.slim_optimize()
        growth_info[genome][predrxn][count]['all_rxns'] = solution

        #turn off all Carbon containing rxns:
        for reaction in model.exchanges:
          for metabolite in reaction.metabolites:
              if contains_carbon(metabolite):
                  # Turn on the reaction by setting a lower bound
                  reaction.lower_bound = 0  # turn off rxn
                  #reaction.upper_bound = 0  # turn off rxn
                  break  # No need to check other metabolites if one already contains carbon

        #save negative control, all rxns off
        solution = model.slim_optimize()
        growth_info[genome][predrxn][count]['no_C'] = solution

        #turn the vitamins back on:
        for vit in vitamins.index:
          try: #turn on given rxn
            model.reactions.get_by_id(vit).lower_bound = -1000
          except KeyError:
            continue

        #save positive control, all vitamins on
        solution = model.slim_optimize()
        growth_info[genome][predrxn][count]['only_vitamins'] = solution

        #turn on the metabolite one by one and then add all single C sources:
        model.reactions.get_by_id(predrxn).lower_bound = -1000
        #optimize with only that compound as a baseline
        solution = model.slim_optimize()
        growth_info[genome][predrxn][count]['no_matti']= solution

        #turn on the first rescue metabolite if there was one, if not save agreement value
        if rescue_comp == 'agreement':
          growth_info[genome][predrxn][count]['first_rescue'] = solution
        else: 
          model.reactions.get_by_id(rescue_comp).lower_bound = -1000
          solution = model.slim_optimize()
          growth_info[genome][predrxn][count]['first_rescue']= solution

        #loop through all the lab metabolites and add back a single C source:
        for exrxn in lab_mets['Exchange_rxn']:
          try: #turn on given rxn
            model.reactions.get_by_id(exrxn).lower_bound = -1000

            #calculate growth rate
            solution = model.slim_optimize()
            growth_info[genome][predrxn][count][exrxn]= solution

            #turn off the rxn again
            model.reactions.get_by_id(exrxn).lower_bound = 0

          except KeyError:
            continue
        
        #turn off the precicted compound again
        model.reactions.get_by_id(predrxn).lower_bound = 0

        #add to count
        count += 1

#Identify percentage of the models and C substrates that produced growth with predicted compound addition

#Prepare some dataframes for storing info
metabolite_recovery = []
metabolite_dict = {}

#We will re-expand the dictionary for each genome
for genome in hq_genomes:
  metabolite_dict[genome] = {}
  for predrxn in growth_info[genome].keys():
    curr_frame = pd.DataFrame.from_dict(growth_info[genome][predrxn],orient="index").T
    sub_vitamins = curr_frame.drop(['all_rxns','no_C','only_vitamins', "no_matti", "first_rescue"]) - curr_frame.loc["first_rescue"].values.squeeze() - curr_frame.loc['no_C'].values.squeeze()
    #print(len(sub_vitamins))

    #Determine whether the models could grow sufficiently (for now we will say a growth rate of >1) with vitamin rescue
    check_growth = sub_vitamins.index[sub_vitamins[sub_vitamins>=1].sum(axis=1)>= (sub_vitamins.shape[1] / 2)].tolist()
    #check_growth = sub_vitamins.index[sub_vitamins[sub_vitamins>=1].any(axis=1)].tolist()
    new_row = [predrxn, genome, len(check_growth),len(check_growth)/len(sub_vitamins)]
    metabolite_recovery.append(new_row)
    metabolite_dict[genome][predrxn] = check_growth
    #print(check_growth)

#metab_recovery_df = pd.DataFrame(metabolite_recovery, columns = ['Predicted_Metabolite', 'Genome','Num_Metabolites','Perc_Metabolites'],index = growth_info[genome].keys())
metab_recovery_df = pd.DataFrame(metabolite_recovery, columns = ['Predicted_Metabolite', 'Genome','Num_Metabolites','Perc_Metabolites'])
#print(metab_recovery_df)
#print(metabolite_dict)
recovered_growth_genomes = [y for y in metab_recovery_df['Num_Metabolites'] if y > 1]
recovered_metabolites = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in metabolite_dict.items()]))
#print(recovered_growth_genomes)
#print(len([y for y in metab_recovery_df['Num_Metabolites'] if y > 0]))

#save the output for number and percentage of genomes containing something:
metab_recovery_df.to_csv("../Output/%s_recovery.csv" %(run_name))
recovered_metabolites.to_csv("../Output/%s_rescued_metabolites.csv" %(run_name))