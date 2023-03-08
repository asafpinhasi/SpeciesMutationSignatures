from ChromosomeMutExtractor import get_mutation_counts
import os
import json
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('scientific_names', help='Specify an input file', type=str)
parser.add_argument('alignment_names', help='Specify an output file', type=str)
args = parser.parse_args()

species_scientific = args.scientific_names.replace('_', ' ').split(";")
species = args.alignment_names.replace('_', ' ').split(";")

chr_file = os.path.join(os.getcwd(), 'chr_files', "chr1.maf.gz")

mutation_counts = get_mutation_counts(chr_file, species)

mutation_folder = 'full_mutation_files'

if not os.path.exists(mutation_folder):
    os.mkdir(mutation_folder)
    
with open(os.path.join(mutation_folder, f"{species_scientific[0]}({species_scientific[1]}).json"),'w') as file:
    json.dump(mutation_counts[0], file)
with open(os.path.join(mutation_folder, f"{species_scientific[1]}({species_scientific[0]}).json"),'w') as file:
    json.dump(mutation_counts[1], file)