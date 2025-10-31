import numpy as np
import pandas as pd
import optparse
import os
import seaborn as sns
import matplotlib.pyplot as plt

#Set up argument parsing
if __name__ == '__main__':
	parser = optparse.OptionParser()
	parser.add_option('-d', '--directory', type='string', dest='file_dir',
					  help="provide the directory that contains the model files")
	parser.add_option('-w', '--work-directory', type='string', dest='work_dir',
					  help="provide the working directory where the compiled output should go")
	#Convert input opts and args to python variables
	opts, args = parser.parse_args()

	file_dir = opts.file_dir
	work_dir = opts.work_dir

	if not file_dir or not work_dir:
		parser.error("Both directory (-d) and work-directory (-w) are required")

#This script will take files from a directory and extract their model consensus
#Then we'll add it to a csv file for analysis
input_dir = os.path.join(file_dir)  # Use file_dir directly as it's already the input path
if not os.path.exists(input_dir):
    raise ValueError(f"Input directory does not exist: {input_dir}")

# Create output directory if it doesn't exist
os.makedirs(work_dir, exist_ok=True)

save_data = [["Genome", "Consensus"]]  # Initialize results list
# Process each file in the input directory
for filename in os.listdir(input_dir):
	if not filename.endswith('_rxn-info.csv'):
		continue
	
	# Get full file path and genome name
	file_path = os.path.join(input_dir, filename)
	genome = filename.replace('_rxn-info.csv', '')
	
	try:
		# Read and process the CSV file
		curr_data = pd.read_csv(file_path)
		per_rxn_consensus = curr_data.sum(numeric_only=True, axis=1) / 60
		avg_consensus = per_rxn_consensus.mean()
		data_row = [genome, avg_consensus]
		save_data.append(data_row)
	except Exception as e:
		print(f"Error processing {filename}: {str(e)}")

# Create DataFrame and save to CSV
save_df = pd.DataFrame(save_data[1:], columns=save_data[0])  # Use first row as column names
output_csv = os.path.join(work_dir, "quality-data.csv")
save_df.to_csv(output_csv, index=False)

# Create and save quality histogram
plt.figure(figsize=(10, 6))
quality_hist = sns.histplot(data=save_df, x="Consensus", kde=True)
plt.title("Model Quality Distribution")
plt.xlabel("Consensus Score")
plt.ylabel("Count")

# Save plot
output_plot = os.path.join(work_dir, "model-quality-histogram.png")
plt.savefig(output_plot)
plt.close()