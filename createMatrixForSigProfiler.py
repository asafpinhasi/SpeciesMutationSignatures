import os
import json
import pandas as pd

def load_dicts_from_folder(folder_path):
    dicts = []
    names = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".json"):
            with open(os.path.join(folder_path, filename), "r") as f:
                dictionary = json.load(f)
                names.append(filename.split('.')[0])
                dicts.append(dictionary)
    return dicts, names

def format_mutation_dict(original_dict):
    new_dict = {}
    for key, value in original_dict.items():
        new_key = key[0] + "[" + key[1] +">" +key[-1] + "]" + key[2]
        new_dict[new_key] = value
    return new_dict


species_mut_count, species_names = load_dicts_from_folder("mutation_files")
new_species_mutation = [format_mutation_dict(species) for species in species_mut_count]
df = pd.DataFrame(new_species_mutation,index=species_names)

file_name = "test_species_mutations.SBS96.txt"
df.T.sort_index().to_csv(file_name, sep = '\t')