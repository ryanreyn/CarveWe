# this version is edited to run on the cluster, it's then coppied and pasted into a text file
# this version is edited to run generally for a single file that has a genome name input as an agument
# run like this: 'python thisscript.py genome_name'

import numpy as np
import pandas as pd
import optparse

from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import pfba
from cobra.medium import minimal_medium
import os 
import sys
import csv
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
  min_components = minimal_medium(model, solution, minimize_components=True)

  #this minimal media function maximizes flux but does not minimize number of components
  max_flux = minimal_medium(model, solution, minimize_components=False)

  #combine min and min1 into a dataframe
  max_flux_df = max_flux.to_frame()
  max_flux_df.rename(columns={0:'max_flux'}, inplace=True)
  min_df = pd.concat([max_flux_df, min_components], axis=1)
  min_df

  return min_df

#this functin adds metabolite names to a dataframe with metabolite ids as the index
def add_met_names(model, min_df):
  #add metabolite names in human readable format
  min_df_names = min_df.copy(deep=True)
  for index in min_df_names.index:
    met = index.replace("EX_", "")
    name = model.metabolites.get_by_id(met).name
    min_df_names.loc[index, "Met_Name"] = name

  return min_df_names

def restore_rxn_states(model, filtered_reactions, rxn_dict, original_bounds):
    for i, rxn_id in enumerate(filtered_reactions):
        rxn = rxn_dict[rxn_id]
        rxn.lower_bound, rxn.upper_bound = original_bounds[i]
  

def test_growth_sensitivity_to_reaction(model, reactions_to_test, verbose=True):
    """
    For each reaction in the provided list, temporarily knock it out and test if growth fails.
    
    Parameters:
        model (cobra.Model): The metabolic model.
        reactions_to_test (list of str): List of reaction IDs to knock out.
        verbose (bool): Whether to print full debug info.
        
    Returns:
        list of str: List of reactions that, when knocked out, cause growth to fail.
    """
    failed_reactions = []

    for rxn_id in reactions_to_test:
        try:
            with model as temp_model:
                if rxn_id not in temp_model.reactions:
                    if verbose:
                        print(f"[SKIP] Reaction {rxn_id} not found in model.")
                    continue

                rxn = temp_model.reactions.get_by_id(rxn_id)
                rxn.knock_out()
                growth = temp_model.slim_optimize()

                if verbose:
                    print(f"[TEST] {rxn_id} knocked out → growth: {growth:.6f}")

                if growth is None or growth <= 1e-6:
                    failed_reactions.append(rxn_id)
                    if verbose:
                        print(f"[FAIL] {rxn_id} knockout caused growth failure.")
        except Exception as e:
            if verbose:
                print(f"[ERROR] Exception testing {rxn_id}: {e}")
            failed_reactions.append(rxn_id)

    return failed_reactions


#Defining a function to convert reaction info csv files to a dictionary
def load_ensemble_map(csv_path):
    #"""Read reaction_id -> ENSEMBLE_STATE from rxn-info.csv file"""
    ensemble_map = {}
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                rxn_id = row[0].strip()
                state = ",".join(cell.strip() for cell in row[1:])  # handle multiple commas
                ensemble_map[rxn_id] = state
    return ensemble_map

def reopen_auto_added_exchanges(model, ensemble_map, lower_bound=-1000., upper_bound=1000., verbose=False):
    """
    Reopen exchange reactions that were auto-added by COBRApy and are not present
    in the ensemble_map. These reactions are assumed to be required for growth.
    
    Parameters:
        model (cobra.Model): The COBRApy model object.
        ensemble_map (dict): Dictionary mapping reaction IDs to ensemble states.
        lower_bound (float): Lower bound to assign to whitelisted reactions.
        upper_bound (float): Upper bound to assign to whitelisted reactions.
        verbose (bool): If True, print each reopened reaction.
    """
    cobra_rxns = set(rxn.id for rxn in model.reactions)
    ensemble_rxns = set(ensemble_map.keys())
    auto_added_rxns = cobra_rxns - ensemble_rxns

    reopened = 0
    for rxn_id in auto_added_rxns:
        if rxn_id.startswith("EX_"):
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = float(lower_bound)
            rxn.upper_bound = float(upper_bound)
            reopened += 1
            if verbose:
                print(f"[AUTO] Reopened auto-added exchange: {rxn.id}")
    
    if verbose:
        print(f"[INFO] Reopened {reopened} auto-added exchange reactions.")

    if verbose:
        print(f"[DEBUG] Final exchange reactions open after whitelisting:")
        for rxn in model.reactions:
            if rxn.id.startswith("EX_") and rxn.lower_bound < 0:
                print(f"  {rxn.id} : lb={rxn.lower_bound}, ub={rxn.upper_bound}")


#We need to import the .xml file or SBML model as well as the reaction states file
model = read_sbml_model('%s/xml_files/%s.xml' %(work_dir, genome))
rxn_info_path = ('%s/ensemble_rxn-info_files/%s_rxn-info.csv' %(work_dir, genome))
rxn_states = pd.read_csv(rxn_info_path, header = None, index_col=0)

#Hard set model solver, we will reinstantiate this on each ensemble run too
model.solver = "cplex"

#Load the ensemble mapping for externally identified reactions to model object reactions
ensemble_map = load_ensemble_map(rxn_info_path)

#Explicitly ensure that the model compartments are correctly defined both in the model
#object and the individual metabolite objects
model.compartments = {
    'c': 'cytosol',
    'p': 'periplasm',
    'e': 'extracellular space'
}

# Also fix all metabolites’ compartment labels
for met in model.metabolites:
    if met.compartment == 'C_c':
        met.compartment = 'c'
    elif met.compartment == 'C_p':
        met.compartment = 'p'
    elif met.compartment == 'C_e':
        met.compartment = 'e'


headers = pd.DataFrame()
media_recipes = pd.DataFrame()

#COBRApy has introduced a new functionality that adds reactions as it reads the SBML
#model file. To contravert this, we will be filtering these out with rxn_info.csv files
filtered_reactions = list(rxn_states.index)

print(model.summary())

print(f"[INFO] Removed {len(model.reactions)-len(filtered_reactions)} reactions from model.")

#Establish a dictionary that corresponds between reaction IDs and COBRA reaction
#objects
rxn_dict = {rxn.id: rxn for rxn in model.reactions}

for j in range(num_models): #column indexing variable
  i = 0 #reaction indexing variable

  if verbose:
    print(f"[INFO] Processing model {j+1}/{num_models}")
    
  # Define the full list of reactions with state=0 in your ensemble file
  zeroed_reactions = [rxn for i, rxn in enumerate(filtered_reactions)
                      if int(rxn_states.iloc[i, j]) == 0]  # test column 0 for simplicity

  # Run the sensitivity test
  problematic = test_growth_sensitivity_to_reaction(
    model, zeroed_reactions, verbose=True
    )

  print("\nReactions that block growth when knocked out:")
  for r in problematic:
      print(f"  - {r}")

  with model as temp_model:

    # Apply knockouts based on ensemble
    for i, rxn_id in enumerate(filtered_reactions):
      state = int(rxn_states.iloc[i, j])
      if rxn_id == "EX_o2_e":
        print(f"[CHECK] EX_o2_e state at ensemble {j}: {state}")
      if state == 0:
          try:
              if rxn_id not in problematic:
                temp_model.reactions.get_by_id(rxn_id).knock_out()
                if verbose:
                    print(f"[KO] Reaction {rxn_id} knocked out.")
          except KeyError:
              if verbose:
                  print(f"[WARN] Reaction {rxn_id} not found in model.")

    #Reopen any auto-added exchange reactions
    reopen_auto_added_exchanges(model, ensemble_map, verbose=False)

    # #Check to see the boundary conditions for exchange reactions
    # print("\n[INFO] Exchange Reaction Bounds")
    # print("=" * 40)
    # for rxn in model.reactions:
    #     if rxn.id.startswith("EX_"):
    #         print(f"{rxn.id:25s} : lb = {rxn.lower_bound:8.2f}, ub = {rxn.upper_bound:8.2f}")

    #Add a couple of lines to try to force dump and reinstantiate solver object
    temp_model.solver = "cplex"

    growth = temp_model.slim_optimize()
    print(f"[ENSEMBLE {j}] Growth rate after knockouts: {growth}")

    #generate minimal media recipe
    try:
      # min_media = find_min_medias(model)
      min_media = find_min_medias(temp_model)

      if min_media.empty:
          print(f"[WARNING] min_media returned empty for model {j}")

      if verbose:
        print(f"[DEBUG] Model {j} returned minimal media with shape: {min_media.shape}")
        print(f"[DEBUG] Non-zero components in minimal media:")
        nonzero_mets = min_media[(min_media != 0).any(axis=1)]
        for met_id in nonzero_mets.index:
            row = nonzero_mets.loc[met_id]
            print(f"  {met_id}: {row.values}")
      
      media_recipes = pd.concat([media_recipes, min_media], axis = 1)

      #generate a header to save the metadata
      model_name = []
      model_number = []
      media_type = []
      ensemble_size = []
      name = genome
      for column in min_media.columns:
        c = 0
        if column == "max_flux":
          c = 1
        media_type.append(c)
        model_name.append(name)
        model_number.append(j)
        ensemble_size.append(num_models)
      header = pd.DataFrame([model_name, model_number, media_type, ensemble_size], columns = min_media.columns, index = ['model_name', 'model_number', 'media_type', 'ensemble_size'])
      #append header info to headers dataframe
      headers = pd.concat([headers, header], axis=1)
    
    except AttributeError:
      print('Media could not be found for %s ensemble %d' %(genome, j))

#add metabolite names
all_media = add_met_names(model, media_recipes)

all_media.to_csv("%s/%s_ensemble_media.csv" %(media_dir, genome))
headers.to_csv("%s/%s_ensemble_header.csv" %(media_dir, genome))
