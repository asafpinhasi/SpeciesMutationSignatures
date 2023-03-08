import pandas as pd
import re
import numpy as np
import gzip
import concurrent.futures
import os
import pickle as pk

from collections import defaultdict
from urllib.request import urlopen

import io

#single chromosome mutation counts extractor
class Chromosome_MutExt: 
    
    
    """
    This class extracts in practice the mutation counts from a given chromosome
    
    Important Attributes:
    chosen_species (String array): scientific name of the 3 species to be used (species[2]) is the outgroup.
    folder (String): The path of the folder containing alignment file
    file_name (String): The file name that the sequencing data is in.
    mutation_count(dict): contains the mutation count extracted from the file
    """
    def __init__(self, 
                 file_name = 'chr1.maf.gz',
                 species = ['oviAri3', 'capHir1', 'bosTau8']):
        
        self.chosen_species = species
        self.file_name = file_name
        self.species1_mutation_count = defaultdict(int)
        self.species2_mutation_count = defaultdict(int)
        
    
    def get_mutation_counts(self):
        
        '''
        Parses the file line by line to avoid reading too much into memory
        Divides the lines into chunks (according to alignment file structure),
        And processes each chunk using __parse_paragraph

        Returns:
            The final count from the whole file
        '''
        
        species = ["s " + s for s in self.chosen_species]
        self.species1_mutation_count = defaultdict(int)
        self.species2_mutation_count = defaultdict(int)

        with gzip.GzipFile(self.file_name, 'rb') as file:
            buffer_size = 10 * 1024 * 1024  # 10MB buffer
            reader = io.BufferedReader(file, buffer_size=buffer_size)
            decoded_reader = io.TextIOWrapper(reader, encoding='utf-8')
            chunk = []
            #n = 10000000
            #i = 0
            for line in decoded_reader:
                #i+=1
                #if i > n:
                #    break
                if line.startswith(tuple(species)):
                    chunk.append(line)
                elif line[0] == 'a': # beggining of new chunk
                    if self.__species_in_chunk(species, chunk):
                        self.__parse_paragraph(chunk, species)
                    chunk = []
        # merges the 192 categories to 96
        self.species1_final_mutation_count = self.__get_final_count(self.species1_mutation_count)
        self.species2_final_mutation_count = self.__get_final_count(self.species2_mutation_count)
        return (self.species1_final_mutation_count, self.species2_final_mutation_count)
     
    
    def __species_in_chunk(self, three_species, chunk):
        for species in three_species:
            if not any(line.startswith(species) for line in chunk):
                return False
        return True
           
    
    anti_pattern = r"[^ATGC]"
    name_regex = re.compile(r"s \w*.")
    sequence_regex = re.compile(r"[atgcnATGCN-]+$")
    
                   
    
    def __parse_paragraph(self, paragraph, species):
        '''
        Parses a single chunk containing sequencing data from the 3 species.
        Extracts from it the mutation count and adds it to mutation_count dict
        '''
        sequences = {}
        for line in paragraph:
            species_name = line[:line.index(".")]
            sequence = self.sequence_regex.findall(line)
            if sequence:
                sequences[species_name] = sequence[0].upper()
            else:
                print("failed:")
                print(line)
                break;
        return self.__count_mutations(sequences, species)

    

    def __count_mutations(self, sequences, species):
        '''
        Counts from the 3 sequences all the mutations that can be detected
        Counts only mutations between the two sister taxas, in which the outlies has
        a nucleotide that is similar to that of one of the taxas.
        Adds the counts to mutation_count
        '''
        
        species1_mutation_count = self.species1_mutation_count
        species2_mutation_count = self.species2_mutation_count
        df = pd.DataFrame(columns = species)
        df[species[0]] = list(sequences[species[0]])
        df[species[1]] = list(sequences[species[1]])
        df[species[2]] = list(sequences[species[2]])
        df[species[0] + " neighbors"] = self.__get_neighbors(sequences[species[0]])
        df[species[1] + " neighbors"] = self.__get_neighbors(sequences[species[1]])
        df[species[2] + " neighbors"] = self.__get_neighbors(sequences[species[2]])
        # creates a boolean column of all point mutation between sister taxa (the nucleotide is different and is not a gap)
        # includes only mutations with adjacent nucleotides for all 3 taxa
        df["mutations"] = ((df[species[0]] != df[species[1]]) & 
                           (df[species[0] + " neighbors"] != '-') & 
                           (df[species[1] + " neighbors"] != '-') & 
                           (df[species[2] + " neighbors"] != '-'))

        # creates a boolean column with only the mutation that happend in one taxa and not in the other
        # this is important for establishing the type of mutation that occurred
        # df[species[0] + " neighbors"].iloc[df["mutations"] & ((df[species[0]] == df[species[2]])]

        mutations = df.loc[df["mutations"] & (df[species[0]] == df[species[2]])]
        for mut in mutations.index:
            mutation = mutations.loc[mut][species[1] + " neighbors"]
            mutation = mutation[0] + mutations.loc[mut][species[2]] + mutation[2] + "->" + mutation[1]
            species1_mutation_count[mutation] += 1

        mutations = df.loc[df["mutations"] & (df[species[1]] == df[species[2]])]
        for mut in mutations.index:
            mutation = mutations.loc[mut][species[0] + " neighbors"]
            mutation = mutation[0] + mutations.loc[mut][species[2]] + mutation[2] + "->" + mutation[1]
            species2_mutation_count[mutation] += 1

        self.species1_mutation_count = species1_mutation_count
        self.species2_mutation_count = species2_mutation_count
        return df

    #Calculated for each position the it's neighboring nucleotides (ignores '-' in alignment)
    def __get_neighbors(self, string):
        nucleotides = {'A','T','G','C'}
        sequence = re.sub(self.anti_pattern, '', string)
        sequence_len = len(sequence)
        i = 0
        neighbors = []
        for idx, char in enumerate(string):
            if i == sequence_len-1 or i == 0:  
                neighbors.append('-')
                if char in nucleotides:
                    i+=1
            elif char not in nucleotides:
                neighbors.append('-')
            else:
                neighbors.append(sequence[i-1:i+2])
                i += 1
        return neighbors

    
    complement = {'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C'}
    
    # merges the mutations properly to reach 96 mutational categories
    def __get_complement(self, mutation):
        comp_mutation = ''.join([self.complement[nuc] for nuc in mutation[2::-1]])
        comp_mutation += "->" + self.complement[mutation[-1]]
        return comp_mutation

    def __get_final_count(self, mutation_count):
        merged_categories = defaultdict(int)

        for mutation, count in mutation_count.items():
            if mutation[1] in {'A', 'G'}:
                merged_categories[self.__get_complement(mutation)] += count
            else:
                merged_categories[mutation] += count
        return merged_categories
    
def get_mutation_counts(file_name, chosen_species):
    chromosome_mutext = Chromosome_MutExt(file_name, chosen_species)
    mutational_counts = chromosome_mutext.get_mutation_counts()
    return mutational_counts
