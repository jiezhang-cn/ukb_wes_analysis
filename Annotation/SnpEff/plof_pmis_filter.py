#!/usr/bin/env python3

import csv
import re
import os
import glob
import multiprocessing as mp
from functools import partial
import time
import operator
from collections import defaultdict

def parse_vcf_header(vcf_file):
    """Parse VCF header to get column names"""
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                return line.strip('#\n').split('\t')
    return None

def extract_ann_fields(info):
    """Extract ANN fields from INFO column"""
    if 'ANN=' not in info:
        return None, None, None
    
    ann_part = info.split('ANN=')[1].split(',')[0]  # Get first ANN annotation
    fields = ann_part.split('|')
    
    if len(fields) < 4:
        return None, None, None
    
    allele = fields[0]
    effect = fields[1]
    gene = fields[3]
    
    return allele, effect, gene

def extract_plof_variants(snpeff_file, output_file):
    """Extract all HIGH impact variants as pLOF variants from SnpEff VCF file"""
    plof_variants = []
    
    with open(snpeff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = fields[1]
            id_field = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Check if this is a HIGH impact variant
            if 'ANN=' in info and '|HIGH|' in info:
                allele, effect, gene = extract_ann_fields(info)
                if allele and effect and gene:
                    # Create formatted variant record
                    variant_id = f"chr{chrom}:{pos}:{ref}:{alt}"
                    plof_variants.append([chrom, pos, variant_id, ref, alt, allele, effect, gene])
    
    # Write pLOF variants to output file in the specified format
    if plof_variants:
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tID\tREF\tALT\tANN[*].ALLELE\tANN[*].EFFECT\tANN[*].GENE\n")
            # Write variant lines
            for variant in plof_variants:
                f.write('\t'.join(variant) + '\n')
    
    return [v[2] for v in plof_variants]  # Return variant IDs

def get_missense_variants(snpeff_file, output_file):
    """Extract all missense variants from SnpEff VCF file"""
    missense_variants = []
    variant_positions = set()  # Use a set for faster lookups
    variant_by_position = {}  # Index variants by position for faster lookup
    
    with open(snpeff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = fields[1]
            id_field = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Check if this is a missense variant
            if 'ANN=' in info and 'missense_variant' in info:
                allele, effect, gene = extract_ann_fields(info)
                if allele and effect and gene:
                    # Create formatted variant record
                    variant_id = f"chr{chrom}:{pos}:{ref}:{alt}"
                    variant_record = [chrom, pos, variant_id, ref, alt, allele, effect, gene]
                    missense_variants.append(variant_record)
                    variant_positions.add(variant_id)
                    
                    # Index by position for faster lookup
                    key = (chrom, pos, ref, alt)
                    variant_by_position[key] = variant_record
    
    # Write missense variants to output file in the specified format
    if missense_variants:
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tID\tREF\tALT\tANN[*].ALLELE\tANN[*].EFFECT\tANN[*].GENE\n")
            # Write variant lines
            for variant in missense_variants:
                f.write('\t'.join(variant) + '\n')
    
    return missense_variants, variant_positions, variant_by_position

def create_position_index(annovar_file):
    """Pre-index the Annovar file by position for faster matching"""
    position_index = defaultdict(list)
    
    with open(annovar_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Read header
        
        # Find relevant column indices
        chr_idx = header.index("Chr")
        pos_idx = header.index("Start")
        ref_idx = header.index("Ref")
        alt_idx = header.index("Alt")
        
        # Read all rows and index by position
        for row_idx, row in enumerate(reader):
            if len(row) <= max(chr_idx, pos_idx, ref_idx, alt_idx):
                continue
            
            chrom = row[chr_idx]
            pos = row[pos_idx]
            
            # Add to position index
            position_index[(chrom, pos)].append(row_idx)
    
    return position_index, header

def filter_deleterious_missense(annovar_file, missense_positions, variant_by_position, output_file):
    """Filter deleterious missense variants predicted by 5+ tools using indexed lookup"""
    deleterious_missense = []
    deleterious_ids = set()  # Use a set for faster lookups
    
    # First, create an index of the Annovar file by position
    print("     Creating position index for Annovar file...")
    position_index_start = time.time()
    position_index, header = create_position_index(annovar_file)
    print(f"     Position index created in {time.time() - position_index_start:.2f}s")
    
    # Find relevant column indices
    chr_idx = header.index("Chr")
    pos_idx = header.index("Start")
    ref_idx = header.index("Ref")
    alt_idx = header.index("Alt")
    sift_idx = header.index("SIFT_pred") if "SIFT_pred" in header else -1
    lrt_idx = header.index("LRT_pred") if "LRT_pred" in header else -1
    pp2_hdiv_idx = header.index("Polyphen2_HDIV_pred") if "Polyphen2_HDIV_pred" in header else -1
    pp2_hvar_idx = header.index("Polyphen2_HVAR_pred") if "Polyphen2_HVAR_pred" in header else -1
    mt_idx = header.index("MutationTaster_pred") if "MutationTaster_pred" in header else -1
    
    # Check if any required column is missing
    if any(idx == -1 for idx in [sift_idx, lrt_idx, pp2_hdiv_idx, pp2_hvar_idx, mt_idx]):
        print("Warning: One or more required prediction columns missing in Annovar file")
        return set(), []
    
    # Now read the Annovar file and only process rows that have matching positions
    with open(annovar_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        
        rows = list(reader)  # Read all rows at once
        
        # Process only variants with matching positions
        for missense_variant in variant_by_position.values():
            chrom = missense_variant[0]
            pos = missense_variant[1]
            ref = missense_variant[3]
            alt = missense_variant[4]
            variant_id = missense_variant[2]
            
            # Skip if we've already processed this variant
            if variant_id not in missense_positions:
                continue
            
            # Look up potential matches in the position index
            position_key = (chrom, pos)
            if position_key not in position_index:
                continue
            
            # Check each potential match
            for row_idx in position_index[position_key]:
                row = rows[row_idx]
                
                # Check if ref and alt match
                if row[ref_idx] != ref or row[alt_idx] != alt:
                    continue
                
                # Check if predicted deleterious by all 5 tools
                sift_pred = row[sift_idx] == "D"
                lrt_pred = row[lrt_idx] == "D"
                pp2_hdiv_pred = row[pp2_hdiv_idx] in ["D", "P"]
                pp2_hvar_pred = row[pp2_hvar_idx] in ["D", "P"]
                mt_pred = row[mt_idx] in ["D", "A"]
                
                if sift_pred and lrt_pred and pp2_hdiv_pred and pp2_hvar_pred and mt_pred:
                    deleterious_missense.append(missense_variant)
                    deleterious_ids.add(variant_id)
                    break  # No need to check other potential matches
    
    # Write deleterious missense variants to output file in the specified format
    if deleterious_missense:
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tID\tREF\tALT\tANN[*].ALLELE\tANN[*].EFFECT\tANN[*].GENE\n")
            # Write variant lines
            for variant in deleterious_missense:
                f.write('\t'.join(variant) + '\n')
    
    return deleterious_ids, deleterious_missense

def merge_variants(plof_file, deleterious_missense, output_file):
    """Merge pLOF and deleterious missense variants into a single file"""
    merged_variants = []
    
    # Read pLOF variants if they exist
    if os.path.exists(plof_file) and os.path.getsize(plof_file) > 0:
        with open(plof_file, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                merged_variants.append(line.strip())
    
    # Add deleterious missense variants
    for variant in deleterious_missense:
        merged_variants.append('\t'.join(variant))
    
    # Write merged variants to output file
    if merged_variants:
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tID\tREF\tALT\tANN[*].ALLELE\tANN[*].EFFECT\tANN[*].GENE\n")
            # Write variants
            for variant in merged_variants:
                f.write(variant + '\n')
    
    return len(merged_variants)

def process_chromosome(chrom, snpeff_dir, annovar_dir, output_dir):
    """Process a single chromosome"""
    # Define file paths with absolute paths
    chrom_str = f"chr{chrom}"
    snpeff_file = os.path.join(snpeff_dir, f"ukb_wes_eur_{chrom_str}_SnpEff.vcf")
    annovar_file = os.path.join(annovar_dir, f"ukb_wes_eur_{chrom_str}.hg38_multianno.csv")
    
    # Output files with absolute paths
    plof_file = os.path.join(output_dir, f"ukb_wes_eur_{chrom_str}_plof.tsv")
    missense_file = os.path.join(output_dir, f"ukb_wes_eur_{chrom_str}_missense.tsv")
    deleterious_missense_file = os.path.join(output_dir, f"ukb_wes_eur_{chrom_str}_deleterious_missense.tsv")
    final_output = os.path.join(output_dir, f"ukb_wes_eur_{chrom_str}_plof_pmis.tsv")
    
    print(f"Processing chromosome {chrom}...")
    print(f"  SnpEff file: {snpeff_file}")
    print(f"  Annovar file: {annovar_file}")
    
    # Check if input files exist
    if not os.path.exists(snpeff_file):
        print(f"  Error: SnpEff file not found: {snpeff_file}")
        return None
    if not os.path.exists(annovar_file):
        print(f"  Error: Annovar file not found: {annovar_file}")
        return None
    
    # Start timer for this chromosome
    start_time = time.time()
    
    print(f"  1. Extracting pLOF variants (all HIGH impact variants)...")
    plof_ids = extract_plof_variants(snpeff_file, plof_file)
    plof_time = time.time()
    print(f"     Found {len(plof_ids)} pLOF variants (took {plof_time - start_time:.2f}s)")
    print(f"     Output: {plof_file}")
    
    print(f"  2. Extracting missense variants...")
    missense_variants, missense_ids, variant_by_position = get_missense_variants(snpeff_file, missense_file)
    missense_time = time.time()
    print(f"     Found {len(missense_ids)} missense variants (took {missense_time - plof_time:.2f}s)")
    print(f"     Output: {missense_file}")
    
    print(f"  3. Filtering for deleterious missense variants...")
    deleterious_ids, deleterious_missense = filter_deleterious_missense(annovar_file, missense_ids, variant_by_position, deleterious_missense_file)
    deleterious_time = time.time()
    print(f"     Found {len(deleterious_ids)} deleterious missense variants (took {deleterious_time - missense_time:.2f}s)")
    print(f"     Output: {deleterious_missense_file}")
    
    print(f"  4. Merging pLOF and deleterious missense variants...")
    total_variants = merge_variants(plof_file, deleterious_missense, final_output)
    merge_time = time.time()
    print(f"     Total variants in final output: {total_variants} (took {merge_time - deleterious_time:.2f}s)")
    print(f"     Final output: {final_output}")
    
    total_time = time.time() - start_time
    print(f"  Complete! Final output written to {final_output}")
    print(f"  Total processing time for chromosome {chrom}: {total_time:.2f} seconds")
    
    return {
        "chrom": chrom,
        "count": total_variants,
        "time": total_time,
        "plof_count": len(plof_ids),
        "missense_count": len(missense_ids),
        "deleterious_count": len(deleterious_ids)
    }

def process_chromosomes_batch(chroms, snpeff_dir, annovar_dir, output_dir):
    """Process a batch of chromosomes in parallel"""
    with mp.Pool(processes=len(chroms)) as pool:
        func = partial(process_chromosome, snpeff_dir=snpeff_dir, annovar_dir=annovar_dir, output_dir=output_dir)
        results = pool.map(func, chroms)
    return [r for r in results if r is not None]

def main():
    # Directories with absolute paths
    snpeff_dir = "/mnt/project/Jzhang_data/exome_data/Annotation/SnpEff/output"
    annovar_dir = "/mnt/project/Jzhang_data/exome_data/Annotation/Annovar/output"
    output_dir = "/opt/notebooks"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Define all chromosomes
    all_chroms = list(range(1, 23))  # 1-22
    
    # Process chromosomes in batches of 5
    batch_size = 5
    batch_results = []
    
    # Main timer for whole process
    start_time = time.time()
    
    for i in range(0, len(all_chroms), batch_size):
        batch = all_chroms[i:i+batch_size]
        print(f"\n{'='*80}")
        print(f"PROCESSING BATCH: Chromosomes {batch}")
        print(f"{'='*80}")
        
        batch_start = time.time()
        results = process_chromosomes_batch(batch, snpeff_dir, annovar_dir, output_dir)
        batch_results.extend(results)
        
        batch_time = time.time() - batch_start
        print(f"\n{'='*80}")
        print(f"BATCH COMPLETE: Chromosomes {batch}")
        print(f"Batch processing time: {batch_time:.2f} seconds")
        print(f"{'='*80}")
    
    # Create a summary file with absolute path
    summary_file = os.path.join(output_dir, "variant_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Chromosome\tTotal Variants\tpLOF Count\tMissense Count\tDeleterious Missense Count\tProcessing Time (s)\n")
        total_all = 0
        total_plof = 0
        total_missense = 0
        total_deleterious = 0
        
        for result in sorted(batch_results, key=lambda x: x["chrom"]):
            f.write(f"chr{result['chrom']}\t{result['count']}\t{result['plof_count']}\t{result['missense_count']}\t{result['deleterious_count']}\t{result['time']:.2f}\n")
            total_all += result['count']
            total_plof += result['plof_count']
            total_missense += result['missense_count']
            total_deleterious += result['deleterious_count']
            
        f.write(f"Total\t{total_all}\t{total_plof}\t{total_missense}\t{total_deleterious}\t{time.time() - start_time:.2f}\n")
    
    total_time = time.time() - start_time
    print(f"\nAll chromosomes processed. Total time: {total_time:.2f} seconds")
    print(f"Summary written to {summary_file}")
    
    # Print overall statistics
    print("\nOverall Statistics:")
    print(f"Total pLOF variants: {total_plof}")
    print(f"Total missense variants: {total_missense}")
    print(f"Total deleterious missense variants: {total_deleterious}")
    print(f"Total variants in final output: {total_all}")

if __name__ == "__main__":
    main()
