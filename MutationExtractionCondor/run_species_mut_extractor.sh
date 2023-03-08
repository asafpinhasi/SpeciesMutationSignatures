#!/bin/bash

# Set the location of the anaconda python
export PATH="/Local/md_keren/anaconda3/bin:$PATH"

# Navigate to the script directory
cd "/storage/md_keren/asafpi/cancer_bioinformatics"

# Run the python script with the specified parameters
/Local/md_keren/anaconda3/bin/python run_chromosome_extraction.py "$@"