#!/usr/bin/env python3
# Changing into the format that SAIGE could recognize

import pandas as pd
import numpy as np
import os
import csv
import re

def extract_gene_from_filename(filename):
    """Extract chromosome number from filename"""
    match = re.search(r'chr(\d+)_plof_pmis\.tsv', filename)
    if match:
        return match.group(1)
    return None

def create_saige_format_file(input_file, output_file):
    """Convert the input file to SAIGE format"""
    print(f"Processing {input_file}...")
    
    try:
        # Read input file
        data = pd.read_csv(input_file, sep='\t')
        
        # Check if necessary columns exist
        if 'ANN[*].GENE' not in data.columns:
            print(f"Error: Column 'ANN[*].GENE' not found in {input_file}")
            return
        
        # Extract all unique genes
        genes = data['ANN[*].GENE'].unique()
        print(f"Found {len(genes)} unique genes")
        
        # Create SAIGE format data
        saige_data = []
        
        for gene in genes:
            # Filter data for current gene
            gene_data = data[data['ANN[*].GENE'] == gene]
            
            # Create 'var' row
            var_row = [gene, 'var']
            var_row.extend(gene_data['ID'].tolist())
            
            # Create 'anno' row
            anno_row = [gene, 'anno']
            
            # Map effects to 'missense' or 'lof'
            # Modified mapping: only missense_variant -> missense, all others -> lof
            effects = []
            for effect in gene_data['ANN[*].EFFECT']:
                if effect == 'missense_variant':
                    effects.append('missense')
                else:
                    effects.append('lof')
            
            anno_row.extend(effects)
            
            # Add rows to SAIGE data
            saige_data.append(var_row)
            saige_data.append(anno_row)
        
        # Write to output file
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            for row in saige_data:
                writer.writerow(row)
        
        print(f"Successfully created {output_file}")
        
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")

def main():
    # Get current directory
    current_dir = os.getcwd()
    
    # Get all files in current directory
    files = os.listdir(current_dir)
    
    # Filter for plof_pmis files
    plof_pmis_files = [f for f in files if f.endswith('_plof_pmis.tsv')]
    
    print(f"Found {len(plof_pmis_files)} plof_pmis files")
    
    # Process each file
    for input_file in plof_pmis_files:
        # Extract chromosome number
        chrom = extract_gene_from_filename(input_file)
        if chrom:
            # Create output filename
            output_file = f"SnpEff_gene_group_chr{chrom}.txt"
            
            # Create SAIGE format file
            create_saige_format_file(os.path.join(current_dir, input_file), 
                                    os.path.join(current_dir, output_file))

if __name__ == "__main__":
    main()
