from Bio import Entrez
from Bio import SeqIO
import os
import re
import datetime
import time
import argparse

def standardize_country(location):
    # Dictionary mapping provinces/regions to their countries
    province_to_country = {
        # Nigeria regions
        'Edo': 'Nigeria',
        'Ondo': 'Nigeria',
        'Nasarawa': 'Nigeria',
        'Ebonyi': 'Nigeria',
        'Lagos': 'Nigeria',
        'Plateau': 'Nigeria',
        
        # Sierra Leone regions
        'Kenema': 'Sierra_Leone',
        'Eastern_Province': 'Sierra_Leone',
        'Western_Area': 'Sierra_Leone',
        
        # Guinea regions
        'Faranah': 'Guinea',
        'Nzerekore': 'Guinea',
        'Gueckedou': 'Guinea',
        
        # Liberia regions
        'Lofa': 'Liberia',
        'Nimba': 'Liberia',
        'Montserrado': 'Liberia'
    }
    
    # Dictionary for standardizing country names
    country_standardization = {
        'Guinea': 'Republic_of_Guinea',
        'Great_Britain': 'United_Kingdom',
        'UK': 'United_Kingdom',
        'England': 'United_Kingdom',
        'Britain': 'United_Kingdom'
    }
    
    # Clean up location string
    location = location.replace(" ", "_")
    location = location.split(',')[0].split(':')[0].strip()
    
    # Check if it's a province we know
    if location in province_to_country:
        return province_to_country[location]
        
    # Check if it's a country name that needs standardization
    if location in country_standardization:
        return country_standardization[location]
    
    return location

def get_reference_length(virus_type):
    # Reference lengths for different Lassa virus segments
    reference_lengths = {
        'L': 7279,  # L segment
        'S': 3402   # S segment
    }
    return reference_lengths.get(virus_type, 7279)  # Default to L segment length

def get_location(record):
    for feature in record.features:
        if feature.type == "source":
            if 'country' in feature.qualifiers:
                return standardize_country(feature.qualifiers['country'][0])
            elif 'geo_loc_name' in feature.qualifiers:
                return standardize_country(feature.qualifiers['geo_loc_name'][0])
    return "Unknown"

def get_collection_date(record):
    for feature in record.features:
        if feature.type == "source":
            if 'collection_date' in feature.qualifiers:
                date_str = feature.qualifiers['collection_date'][0]
                
                if date_str in ['', 'unknown', 'Unknown', 'NA', 'missing']:
                    return 9999
                
                try:
                    date_str = date_str.replace('_', '-').replace('/', '-')
                    
                    if '-' in date_str:
                        parts = date_str.split('-')
                        potential_year = max(int(parts[0]), int(parts[-1]))
                        if 1900 <= potential_year <= 2024:
                            return potential_year
                    else:
                        year = int(date_str)
                        if 1900 <= year <= 2024:
                            return year
                        
                except (ValueError, IndexError):
                    year_match = re.search(r'(19|20)\d{2}', date_str)
                    if year_match:
                        return int(year_match.group())
                    
    return 9999

def download_sequences(segment_choice, genome_choice, host_choice, metadata_choice, completeness_threshold=0):
    """Download sequences from GenBank"""
    Entrez.email = "anonymous@example.com"
    
    # Reference sequence IDs for each segment
    reference_sequences = {
        'L_segment': 'NC_004297',  # Lassa virus L segment
        'S_segment': 'NC_004296'   # Lassa virus S segment
    }

    # Build query with all filters
    query_parts = []
    
    # Add virus and segment filter
    segment_term = "L segment" if segment_choice == "L_segment" else "S segment"
    query_parts.append(f'("Lassa virus"[organism] OR "Lassa mammarenavirus"[organism]) AND "{segment_term}"[Title]')

    # Add host filter
    if host_choice == '1':  # Human
        query_parts.append('("Homo sapiens"[host] OR "human"[host]) NOT ("Mastomys"[host] OR "rodent"[host] OR "Unknown"[host])')
    elif host_choice == '2':  # Non-human
        query_parts.append('("Mastomys"[host] OR "rodent"[host]) NOT ("Homo sapiens"[host] OR "human"[host])')
    
    # Add genome completeness filter
    if genome_choice == '1':  # Complete genomes only
        query_parts.append('("complete genome"[Title] OR "complete sequence"[Title])')
    elif genome_choice == '2' and completeness_threshold > 0:  # Partial genomes with threshold
        ref_length = get_reference_length(segment_choice.split('_')[0])
        min_length = int(ref_length * (completeness_threshold/100))
        query_parts.append(f'"{min_length}"[SLEN]:"{ref_length}"[SLEN]')
    
    # Combine all query parts
    query = " AND ".join(query_parts)
    
    # Download sequences
    output_file = "downloaded_genomes.gb"
    ref_id = reference_sequences[segment_choice]
    
    # First, get the list of IDs from the search
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    
    if ref_id in id_list:
        id_list.remove(ref_id)
    
    print(f"\nFound {len(id_list)} sequences matching criteria...")
    
    with open(output_file, "w") as out_handle:
        # Download reference sequence first
        print(f"\nDownloading reference sequence {ref_id}...")
        ref_handle = Entrez.efetch(db="nucleotide", id=ref_id, rettype="gb", retmode="text")
        out_handle.write(ref_handle.read())
        ref_handle.close()
        time.sleep(1)
        
        # Download remaining sequences in batches
        batch_size = 100
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i:i + batch_size]
            print(f"Downloading batch {(i//batch_size)+1} of {(len(id_list)-1)//batch_size + 1}...")
            fetch_handle = Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text")
            out_handle.write(fetch_handle.read())
            fetch_handle.close()
            time.sleep(1)
    
    return output_file, query

def cli_main():
    """Entry point for command line interface"""
    parser = argparse.ArgumentParser(description='Download and analyze Lassa virus sequences')
    
    # Required arguments
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory for results')
    parser.add_argument('--segment', type=str, required=True, choices=['L', 'S'],
                       help='Viral segment to analyze (L or S)')
    
    # Optional arguments
    parser.add_argument('--genome', type=str, choices=['1', '2', '3'],
                       help='Genome type: 1=Complete, 2=Partial, 3=All')
    parser.add_argument('--completeness', type=float,
                       help='Minimum completeness percentage (1-100) when using --genome 2')
    parser.add_argument('--host', type=str, choices=['1', '2', '3'],
                       help='Host: 1=Human, 2=Non-human, 3=All')
    parser.add_argument('--metadata', type=str, choices=['1', '2', '3', '4'],
                       help='Metadata filter: 1=Location, 2=Date, 3=Both, 4=None')
    parser.add_argument('--beast', type=str, choices=['1', '2'],
                       help='BEAST format: 1=No, 2=Yes')
    
    args = parser.parse_args()
    
    # Convert paths to absolute paths
    args.output_dir = os.path.abspath(args.output_dir)
    
    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(args.output_dir)
    
    try:
        # Convert segment choice to internal format
        segment_choice = f"{args.segment}_segment"
        
        # Check if all required parameters are provided
        non_interactive = all([args.genome, args.host, args.metadata])
        if args.genome == '2' and args.completeness is None:
            non_interactive = False
        if args.metadata in ['2', '3'] and args.beast is None:
            non_interactive = False
            
        # Run pipeline
        download_sequences(segment_choice, args.genome, args.host, args.metadata, args.completeness)
        
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    cli_main() 