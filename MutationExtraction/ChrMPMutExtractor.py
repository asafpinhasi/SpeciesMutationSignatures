#multiprocessing over one species triplet and many chromosomes
import pandas as pd
import re
import numpy as np
import gzip
import concurrent.futures
import os
import pickle as pk
import sys
import json

import multiprocessing

from collections import defaultdict
from urllib.request import urlopen
from ChromosomeMutExtractor import get_mutation_counts
'''
To run extraction on your computer you first need to download the files from:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/.

Run the commads:
mut_extractor = MPMutExtractor(folder, file_names, species, url)
mut_extractor.get_mutation_counts()

Where:
folder: The folder containing the downloaded maf files with the different chromosomes alignment data
file_names: The name of all the files with the data from the folder
species: A String array with the three (scientific) names of the species you want to analyse.
example = ["Ovis aries", "Capra hircus", "Bos taurus"], where Bos taurus is the outgroup
url - no need to specify, but should contain the url of the multiz100way in UCSC
'''



# Multiprocessing mutational count extractor
class MPMutExtractor:
    """
    This class uses multiprocessing to extract mutation counts
    It extracts counts from multiple chromosome files concurrently using Chromosome_MutSig class.
    
    Important Attributes:
    species_names (String array): scientific name of the 3 species to be used (species[2]) is the outlgroup.
    folder (String): The path of the folder containing alignment files
    file_names (String array): The file names that the sequencing files are in.
    """
    

    def __init__(self, 
                 folder = r"C:\Users\KerenYlab\Downloads" + "\\", 
                 chromosomes_folder = "chr_files",
                 species = ["Ovis aries", "Capra hircus", "Bos taurus"],
                 url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/"):
        
        """
        A constructor to initialize a new MPMutExtractor instance.
        
        Parameters:
        self.html : html page information with about maf alignment files
        folder, file_names, species_names : as described above
        """
        
        
        page = urlopen(url)
        html_bytes = page.read()
        self.html = html_bytes.decode("utf-8")
        self.folder = os.path.join(os.getcwd(), chromosomes_folder)
        self.chromosomes = [file.split('.')[0] for file in os.listdir(chromosomes_folder) if file.startswith("chr")]
        chromosome_files = [os.path.join(self.folder, file) for file in os.listdir(chromosomes_folder) if file.startswith("chr")]
        self.file_names = chromosome_files
        self.species_names = species
    
    
    def get_mutation_counts(self, output_dir = "mutation_files"):
        """
        This function that runs the concurrent mutation count.

        Returns:
        a list of the mutation count for each chromosome (each mutation count is a dict)
        """
        self.get_species_alignment_names(self.species_names)     
        results = None
        # Create an Executor object that runs the processes
        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
            results = executor.map(get_mutation_counts, [file for file in self.file_names], [self.chosen_species for i in self.file_names])
            self.mutation_counts = list(results)
        
        species = self.species_names
        mutation_counts = self.mutation_counts
        chromosomes = self.chromosomes

        
        output_dir1 = f"{output_dir}_{species[0]}({species[1]})"
        output_dir2 = f"{output_dir}_{species[1]}({species[0]})"

        if not os.path.exists(output_dir1):
            os.mkdir(output_dir1)

        if not os.path.exists(output_dir2):
            os.mkdir(output_dir2)

        for i in range(len(mut_extractor.chromosomes)):
            with open(os.path.join(output_dir1, f"{chromosomes[i]}.json"),'w') as file:
                json.dump(mutation_counts[i][0], file)
            with open(os.path.join(output_dir2, f"{chromosomes[i]}.json"),'w') as file:
                json.dump(mutation_counts[i][1], file)
        
        return self.mutation_counts

            
        
    # finds scientific names
    def __find_scientific_names(self):
        start = re.search("Assemblies used in these alignments:", self.html)
        end = re.search("---------------------------------------------------------------", self.html)
        text = self.html[start.end():end.start()]
        names = []
        for line in text.splitlines():
            if (line != '') and not line.isspace() and line[0] != '=':
                names.append(line[25:56].strip(' '))
        self.scientific_names = names
        return names


    #finds names for alignment
    def __find_alignment_names(self):
        names = []
        start = re.search("Assemblies used in these alignments:", self.html)
        end = re.search("---------------------------------------------------------------", self.html)
        text = self.html[start.end():end.start()]
        for line in text.splitlines():
            if (line == '') or line.isspace() or line[0] == '=':
                continue
            names.append(re.findall("/\S*\s", line)[-1][1:-1])
        self.alignment_names = names
        return names
    
    #find common names
    def __find_common_names(self):
        start = re.search("Assemblies used in these alignments:", self.html)
        end = re.search("---------------------------------------------------------------", self.html)
        text = self.html[start.end():end.start()]
        names = []
        for line in text.splitlines():
            if (line != '') and not line.isspace() and line[0] != '=':
                names.append(line[:25].strip(' '))
        self.common_names = names
        return names
    
    #gets the alignment names for the 3 species given to it
    def get_species_alignment_names(self, species = ["Ovis aries", "Capra hircus", "Bos taurus"]):
        
        if not hasattr(self, "scientific_names") or not hasattr(self, "alignment_names"):
            self.common_names = self.__find_common_names()
            self.scientific_names = self.__find_scientific_names()
            self.alignment_names = self.__find_alignment_names()
        
        scientific_names = self.scientific_names
        alignment_names = self.alignment_names
        common_names = self.common_names
        
        species = self.species_names
        chosen_species = [None]*3

        # 3 chosen taxa = Ovis aries, Capra hircus and outgroup Bos taurus
        chosen_species[0] = np.array(alignment_names)[(np.array(scientific_names) == species[0])][0]
        chosen_species[1] = np.array(alignment_names)[(np.array(scientific_names) == species[1])][0]
        chosen_species[2] = np.array(alignment_names)[(np.array(scientific_names) == species[2])][0]

        self.chosen_species = chosen_species
        return chosen_species
    

# species = ['Homo sapiens', 'Pan troglodytes', 'Gorilla gorilla gorilla']

# mut_extractor = MPMutExtractor(species = species)
# mut_extractor.get_species_alignment_names()
# mutation_count = mut_extractor.get_mutation_counts('mutation_files_test')