from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm
import argparse
import os
import re
import shutil
import tempfile
from datetime import datetime
import calendar
import time
import sys
from io import StringIO
from collections import defaultdict, OrderedDict
import subprocess
import site

# Add these constants at the top of the file
REFERENCE_SEQUENCES = {
    'L': {
        'id': 'NC_004297.1',
        'description': 'Lassa virus segment L, complete sequence',
        'location': 'SierraLeone',
        'date': 'NA'
    },
    'S': {
        'id': 'NC_004296.1',
        'description': 'Lassa virus segment S, complete sequence',
        'location': 'SierraLeone',
        'date': 'NA'
    }
}

PINNEO_SEQUENCES = {
    'L': {
        'id': 'KM822127.1',
        'description': 'Lassa virus strain Pinneo segment L, complete sequence',
        'location': 'Nigeria',
        'date': '1969.000'
    },
    'S': {
        'id': 'KM822128.1',
        'description': 'Lassa virus strain Pinneo segment S, complete sequence',
        'location': 'Nigeria',
        'date': '1969.000'
    }
}

# Josiah RefSeq CDS coordinates (1-based inclusive). L and GPC are complement genes.
# NC_004297.1 (L segment): Z 66..365 (+); L complement(466..7122)
# NC_004296.1 (S segment): NP 101..1810 (+); GPC complement(1872..3347)
REFERENCE_COORDS = {
    'L': {
        'id': 'NC_004297.1',
        'L_protein': {'start': 466, 'end': 7122, 'strand': '-'},
        'Z_protein': {'start': 66, 'end': 365, 'strand': '+'},
    },
    'S': {
        'id': 'NC_004296.1',
        'NP': {'start': 101, 'end': 1810, 'strand': '+'},
        'GPC': {'start': 1872, 'end': 3347, 'strand': '-'},
    }
}

# Protein (CDS) alignment options — names match CLI / interactive choices
PROTEIN_NUM = {'1': 'L', '2': 'Z', '3': 'NP', '4': 'GPC'}
PROTEIN_DESCRIPTIONS = {
    'L': 'RNA-dependent RNA polymerase (L segment)',
    'Z': 'RING finger / Z matrix protein (L segment)',
    'NP': 'nucleoprotein (S segment)',
    'GPC': 'glycoprotein precursor (S segment)',
}
LASSA_PROTEIN_REGIONS = {
    'L':   {'segment': 'L', 'ref_id': 'NC_004297.1', 'start': 466, 'end': 7122, 'strand': '-'},
    'Z':   {'segment': 'L', 'ref_id': 'NC_004297.1', 'start': 66, 'end': 365, 'strand': '+'},
    'NP':  {'segment': 'S', 'ref_id': 'NC_004296.1', 'start': 101, 'end': 1810, 'strand': '+'},
    'GPC': {'segment': 'S', 'ref_id': 'NC_004296.1', 'start': 1872, 'end': 3347, 'strand': '-'},
}

def get_segment_type(record):
    """Determine if sequence is L or S segment based on proteins in header"""
    description = record.description.lower()
    
    # First check for direct segment mentions
    if 'segment l' in description or 'l segment' in description:
        return 'L'
    if 'segment s' in description or 's segment' in description:
        return 'S'
    
    # S segment proteins - expanded list
    s_proteins = [
        'nucleoprotein', 
        'glycoprotein',
        'nucleocapsid',
        'gpc',
        'np',
        'gc',
        'gn',
        'g1',
        'g2'
    ]
    
    # L segment proteins - expanded list
    l_proteins = [
        'polymerase',
        'rna-dependent rna polymerase',
        'rna dependent rna polymerase',
        'rdrp',
        'z protein',
        'zinc finger protein',
        'zinc-finger protein',
        'matrix protein',
        'l protein',
        'large protein'
    ]
    
    for protein in s_proteins:
        if protein in description:
            return 'S'
    for protein in l_proteins:
        if protein in description:
            return 'L'
    
    # If no protein found in header, try original header stored during download
    if hasattr(record, 'original_header'):
        description = record.original_header.lower()
        # Check for direct segment mentions in original header
        if 'segment l' in description or 'l segment' in description:
            return 'L'
        if 'segment s' in description or 's segment' in description:
            return 'S'
        
        for protein in s_proteins:
            if protein in description:
                return 'S'
        for protein in l_proteins:
            if protein in description:
                return 'L'
    
    return None

def convert_date_to_decimal_year(date_str):
    """Convert date to decimal year format"""
    try:
        # First, clean the date string
        date_str = date_str.strip()
        
        # If it's just a year
        if date_str.isdigit() and len(date_str) == 4:
            return f"{int(date_str)}.000"
        
        # Try different date formats
        formats = [
            ('%Y-%m-%d', True),    # 2013-08-15 -> full precision
            ('%Y-%m', False),      # 2013-08 -> month precision
            ('%d-%b-%Y', True),    # 15-Aug-2013 -> full precision
            ('%b-%Y', False),      # Aug-2013 -> month precision
        ]
        
        for date_format, is_full_date in formats:
            try:
                date = datetime.strptime(date_str, date_format)
                year = date.year
                
                if is_full_date:
                    # Calculate exact day for full dates
                    day_of_year = date.timetuple().tm_yday
                    total_days = 366 if calendar.isleap(year) else 365
                    decimal_year = year + (day_of_year / total_days)
                    return f"{decimal_year:.3f}"
                else:  # Month precision
                    # Use middle of the month
                    middle_date = datetime(year, date.month, 15)
                    day_of_year = middle_date.timetuple().tm_yday
                    total_days = 366 if calendar.isleap(year) else 365
                    decimal_year = year + (day_of_year / total_days)
                    return f"{decimal_year:.3f}"
            except ValueError:
                continue
        
        return "UnknownDate"
            
    except Exception as e:
        print(f"Date conversion error for {date_str}: {str(e)}")
        return "UnknownDate"

def get_standardized_country_name(country):
    """Convert country name to standardized format"""
    # Dictionary mapping various country name formats to standardized names
    country_mappings = {
        # Sierra Leone variations
        'sierra leone': 'SierraLeone',
        'sierra-leone': 'SierraLeone',
        'sierraleone': 'SierraLeone',
        'Sierra Leone': 'SierraLeone',
        'Sierra-Leone': 'SierraLeone',
        'SIERRA LEONE': 'SierraLeone',
        
        # Nigeria variations
        'nigeria': 'Nigeria',
        'Nigeria': 'Nigeria',
        'NIGERIA': 'Nigeria',
        
        # Liberia variations
        'liberia': 'Liberia',
        'Liberia': 'Liberia',
        'LIBERIA': 'Liberia',
        
        # Guinea variations
        'guinea': 'Guinea',
        'Guinea': 'Guinea',
        'GUINEA': 'Guinea',
        
        # Mali variations
        'mali': 'Mali',
        'Mali': 'Mali',
        'MALI': 'Mali',
        
        # Ghana variations
        'ghana': 'Ghana',
        'Ghana': 'Ghana',
        'GHANA': 'Ghana',
        
        # Benin variations
        'benin': 'Benin',
        'Benin': 'Benin',
        'BENIN': 'Benin',
        
        # Burkina Faso variations
        'burkina faso': 'BurkinaFaso',
        'Burkina Faso': 'BurkinaFaso',
        'BURKINA FASO': 'BurkinaFaso',
        'burkinafaso': 'BurkinaFaso',
        'burkina-faso': 'BurkinaFaso',
        'Burkina-Faso': 'BurkinaFaso',
        
        # Ivory Coast variations
        "cote d'ivoire": 'IvoryCoast',
        "Cote d'Ivoire": 'IvoryCoast',
        "COTE D'IVOIRE": 'IvoryCoast',
        "côte d'ivoire": 'IvoryCoast',
        "Côte d'Ivoire": 'IvoryCoast',
        'ivory coast': 'IvoryCoast',
        'Ivory Coast': 'IvoryCoast',
        'IVORY COAST': 'IvoryCoast',
        'ivorycoast': 'IvoryCoast',
        'ivory-coast': 'IvoryCoast',
        'Ivory-Coast': 'IvoryCoast',
        
        # Togo variations
        'togo': 'Togo',
        'Togo': 'Togo',
        'TOGO': 'Togo'
    }
    
    # First try exact match
    if country in country_mappings:
        return country_mappings[country]
    
    # Then try lowercase match
    return country_mappings.get(country.lower(), country)

def standardize_city_name(city):
    """Standardize city/province names according to rules"""
    if not city or city == "UnknownCity":
        return "UnknownCity"
        
    # Remove any text after comma
    city = city.split(',')[0].strip()
    
    # Special case replacements
    city_replacements = {
        "KGH": "Kenema",
        "Kgh": "Kenema",
        "kgh": "Kenema",
        "N'Zerekore": "Nzerekore",
        "N'Zérékoré": "Nzerekore",
        "N'zerekore": "Nzerekore"
    }
    
    for old, new in city_replacements.items():
        if city == old:
            return new
    
    # Remove common words to exclude
    words_to_remove = [
        "State", "Hospital", "Medical Center", 
        "Medical Centre", "Clinic", "Health Center"
    ]
    for word in words_to_remove:
        city = city.replace(word, "").strip()
    
    # Handle special cases of all caps
    if city.isupper():
        city = city.title()
    
    # Remove double spaces
    city = " ".join(city.split())
    
    return city.strip()

def get_metadata(record):
    """Extract metadata from record"""
    # Get accession - ensure no extra prefixes
    accession = record.id.replace('_R_', '')  # Remove _R_ prefix if present
    
    # Initialize location, province, date, and host
    location = "UnknownLoc"
    province = "UnknownCity"
    collection_date = "UnknownDate"
    host_type = "UnknownHost"
    
    # Look for location, date, and host in source features
    for feature in record.features:
        if feature.type == "source":
            # Get location and province from geo_loc_name
            if 'geo_loc_name' in feature.qualifiers:
                geo_loc = feature.qualifiers['geo_loc_name'][0]
                # Only use location if it's not 'missing'
                if 'missing' not in geo_loc.lower():
                    # Split location into country and province
                    parts = [part.strip() for part in geo_loc.split(':')]
                    location = get_standardized_country_name(parts[0])
                    # Get province if available
                    if len(parts) > 1:
                        province = standardize_city_name(parts[1])
            
            # Get host information
            if 'host' in feature.qualifiers:
                host = feature.qualifiers['host'][0].lower()
                if is_human_host(host):
                    host_type = "Human"
                elif is_rodent_host(host):
                    host_type = "Rodent"
                else:
                    host_type = "Other"
            
            # Get collection date
            if 'collection_date' in feature.qualifiers:
                date_str = feature.qualifiers['collection_date'][0]
                if date_str.lower() != 'missing':
                    collection_date = convert_date_to_decimal_year(date_str)
    
    # Format the metadata string
    if location == "UnknownLoc":
        return f"{accession}_{location}_{host_type}_{collection_date}"
    else:
        return f"{accession}_{location}_{province}_{host_type}_{collection_date}"

def fetch_sequences():
    Entrez.email = "anonymous@example.com"
    max_retries = 3
    retry_delay = 5  # seconds
    
    print("Searching for Lassa virus sequences...")
    search_term = "Mammarenavirus lassaense[Organism]"
    
    # Try the initial search with retries
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            record = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            if attempt == max_retries - 1:  # Last attempt
                print(f"Error: Unable to connect to NCBI after {max_retries} attempts.")
                print(f"Please check your internet connection or try again later.")
                print(f"Error details: {str(e)}")
                raise
            print(f"Connection attempt {attempt + 1} failed. Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
    
    count = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    
    print(f"Found {count} total sequences to process")
    
    batch_size = 100
    sequences = []
    segment_counts = {'L': 0, 'S': 0, 'Unknown': 0}
    location_counts = {'L': {}, 'S': {}}
    
    print("Downloading sequences in batches...")
    for start in tqdm(range(0, count, batch_size), desc="Downloading batches"):
        for attempt in range(max_retries):
            try:
                handle = Entrez.efetch(db="nucleotide", 
                                     rettype="gb", 
                                     retmode="text",
                                     retstart=start,
                                     retmax=batch_size,
                                     webenv=webenv,
                                     query_key=query_key)
                                     
                batch_records = SeqIO.parse(handle, "genbank")
                
                for record in batch_records:
                    # Store original header
                    original_header = record.description
                    record.original_header = original_header
                    
                    segment_type = get_segment_type(record)
                    sequences.append({
                        'id': record.id,
                        'record': record,
                        'segment': segment_type,
                        'original_header': original_header
                    })
                    
                    # Update counts
                    if segment_type:
                        segment_counts[segment_type] += 1
                    else:
                        segment_counts['Unknown'] += 1
                
                handle.close()
                break  # Success, exit retry loop
                
            except Exception as e:
                if attempt == max_retries - 1:  # Last attempt
                    print(f"\nError downloading batch starting at {start}: {str(e)}")
                    print("Skipping this batch and continuing...")
                    continue
                print(f"\nRetry {attempt + 1} for batch starting at {start}...")
                time.sleep(retry_delay)
    return sequences, segment_counts, location_counts

def get_segment_from_protein(description):
    """Helper function to determine segment type from protein names in header"""
    description = description.lower()
    
    # S segment proteins
    s_proteins = ['nucleoprotein', 'glycoprotein', 'gpc', 'np']
    # L segment proteins
    l_proteins = ['polymerase', 'z protein', 'l protein']
    
    for protein in s_proteins:
        if protein in description:
            return 'S'
    for protein in l_proteins:
        if protein in description:
            return 'L'
    
    return None

def process_sequences(sequences, segment_type):
    """Process sequences based on segment type (L, S, or both)"""
    if segment_type.upper() == 'BOTH':
        l_sequences = []
        s_sequences = []
        unknown_sequences = []
        unknown_headers = []
        
        for seq in sequences:
            record = seq['record']
            original_header = seq['original_header']
            
            # First try the official segment annotation
            if seq['segment'] in ['L', 'S']:
                final_segment = seq['segment']
            else:
                # For unknown segments, try to determine from protein names
                final_segment = get_segment_from_protein(original_header)
            
            # Only modify header after classification is complete
            record.id = get_metadata(record)
            record.description = ''
            
            if final_segment == 'L':
                l_sequences.append(record)
            elif final_segment == 'S':
                s_sequences.append(record)
            else:
                unknown_sequences.append(record)
                unknown_headers.append(original_header)
                
        return l_sequences, s_sequences, unknown_sequences, unknown_headers
    else:
        filtered_sequences = []
        target_segment = segment_type.upper()
        
        for seq in sequences:
            record = seq['record']
            original_header = seq['original_header']
            
            # First try the official segment annotation
            if seq['segment'] == target_segment:
                final_segment = seq['segment']
            else:
                # For unknown segments, try to determine from protein names
                final_segment = get_segment_from_protein(original_header)
            
            if final_segment == target_segment:
                record.id = get_metadata(record)
                record.description = ''
                filtered_sequences.append(record)
                    
        return filtered_sequences

def calculate_location_counts(sequences):
    """Calculate location counts for each segment type"""
    locations = {'L': {}, 'S': {}, 'unknown': {}}
    
    for seq in sequences:
        record = seq['record']
        segment = get_segment_type(record)
        segment_type = 'unknown' if segment is None else segment
        
        location = "UnknownLoc"
        for feature in record.features:
            if feature.type == "source":
                if 'geo_loc_name' in feature.qualifiers:
                    geo_loc = feature.qualifiers['geo_loc_name'][0]
                    if 'missing' not in geo_loc.lower():
                        location = geo_loc.split(':')[0].strip()
                        location = get_standardized_country_name(location)
        
        if location not in locations[segment_type]:
            locations[segment_type][location] = 0
        locations[segment_type][location] += 1
    
    return locations

def write_summary(outdir, initial_total, filtered_total, initial_segment_counts, location_counts, requested_segment, written_counts, filtered_sequences, genome_choice, completeness=None, sequences=None, host_choice=None, completeness_filtered=None):
    """Write summary report to a file"""
    summary_file = os.path.join(outdir, "summary_Lassa.txt")
    
    with open(summary_file, 'w') as f:
        f.write("Lassa Virus Sequence Download Summary\n")
        f.write("====================================\n\n")
        
        # Document filtering steps
        write_filtering_steps(f, sequences, genome_choice, completeness, 
                            host_choice, completeness_filtered, filtered_sequences)
        
        # Remove Output Files section from here - it will be written after metadata filtering

def write_geographical_distribution(f, locations, title):
    """Helper function to write geographical distribution tables"""
    f.write(f"\n{title}:\n")
    f.write("-" * len(title) + "\n")
    f.write("Country        L segments    S segments    Unknown         Total\n")
    f.write("-----------   -----------   -----------   ---------   ----------\n")
    
    all_countries = set()
    for segment_type in ['L', 'S', 'unknown']:
        all_countries.update(locations[segment_type].keys())
    
    total_l = 0
    total_s = 0
    total_unknown = 0
    for country in sorted(all_countries):
        l_count = locations['L'].get(country, 0)
        s_count = locations['S'].get(country, 0)
        unknown_count = locations['unknown'].get(country, 0)
        country_total = l_count + s_count + unknown_count
        
        if country_total > 0:
            f.write(f"{country:<13} {l_count:>11}   {s_count:>11}   {unknown_count:>9}   {country_total:>10}\n")
        
        total_l += l_count
        total_s += s_count
        total_unknown += unknown_count
    
    f.write("-" * 60 + "\n")
    total_all = total_l + total_s + total_unknown
    f.write(f"{'Total':<13} {total_l:>11}   {total_s:>11}   {total_unknown:>9}   {total_all:>10}\n\n")

def get_user_input(prompt, valid_options):
    """Helper function to get valid user input"""
    while True:
        print(prompt)
        choice = input("Enter your choice: ").strip()
        if choice in valid_options:
            return choice
        print(f"Invalid choice. Please choose from {', '.join(valid_options)}")

def get_completeness():
    """Get valid completeness percentage"""
    while True:
        try:
            completeness = float(input("Enter minimum sequence completeness (1-100): "))
            if 1 <= completeness <= 100:
                return completeness
            print("Please enter a value between 1 and 100")
        except ValueError:
            print("Please enter a valid number")

def is_complete_sequence(record):
    """Check if sequence is complete based on length (>99% of reference)"""
    reference_lengths = {
        'L': 7279,  # Reference length for L segment
        'S': 3402   # Reference length for S segment
    }
    
    # Get segment type
    segment = get_segment_type(record)
    if not segment:
        return False
        
    # Calculate completeness percentage
    seq_length = len(record.seq)
    ref_length = reference_lengths[segment]
    completeness = (seq_length / ref_length) * 100
    
    # Consider sequence complete if it's at least 99% of the reference length
    return completeness >= 99

def meets_minimum_completeness(record, min_completeness):
    """Check if sequence meets minimum completeness percentage"""
    reference_lengths = {
        'L': 7279,  # Reference length for L segment
        'S': 3402   # Reference length for S segment
    }
    
    # Get segment type
    segment = get_segment_type(record)
    if not segment:
        return False
    
    # Calculate completeness
    seq_length = len(record.seq)
    ref_length = reference_lengths[segment]
    completeness = (seq_length / ref_length) * 100
    
    return completeness >= min_completeness

def calculate_segment_counts(sequences):
    """Calculate segment counts from a list of sequences"""
    counts = {'L': 0, 'S': 0, 'unknown': 0}
    for seq in sequences:
        segment = get_segment_type(seq['record'])
        if segment:
            counts[segment] += 1
        else:
            counts['unknown'] += 1
    return counts

def analyze_hosts(sequences):
    """Print all unique hosts found in sequences"""
    host_counts = {}
    no_host_count = 0
    
    for seq in sequences:
        record = seq['record']
        host_found = False
        
        for feature in record.features:
            if feature.type == "source":
                if 'host' in feature.qualifiers:
                    host = feature.qualifiers['host'][0].lower()
                    if host not in host_counts:
                        host_counts[host] = 0
                    host_counts[host] += 1
                    host_found = True
                    break
        
        if not host_found:
            no_host_count += 1
    
    print("\nHost distribution in sequences:")
    print("-----------------------------")
    for host, count in sorted(host_counts.items(), key=lambda x: (-x[1], x[0])):
        print(f"{host}: {count} sequences")
    if no_host_count > 0:
        print(f"\nSequences with no host information: {no_host_count}")

def is_human_host(host):
    """Check if the host is human"""
    # Collapse spaces so "H. sapiens" and "H.sapiens" both match
    h = host.lower().strip().replace(' ', '')
    human_patterns = [
        'homosapiens',
        'homonsapiens',  # Common misspelling
        'homosapien',    # Common misspelling
        'h.sapiens',     # Common GenBank form
        'hsapiens',
        'human',
        'humans',
        'humanpatient',
        'humanisolate',
        'patient',
    ]
    return any(pattern in h for pattern in human_patterns)

def is_rodent_host(host):
    """Check if the host is a rodent"""
    rodent_patterns = [
        # Specific species
        'mastomys natalensis',
        'mastomys',
        'natal multimammate',
        'multimammate rat',
        'multimammate mouse',
        'hylomyscus pamfi',
        'cavia porcellus',
        'guinea pig',
        # Genera
        'mus',
        'rattus',
        'lophuromys',
        'praomys',
        'hylomyscus',  # Added genus
        'cavia',
        # Common names
        'rodent',
        'mouse',
        'mice',
        'rat',
        # Scientific terms
        'natalensis',
        'murinae',
        'muridae'
    ]
    return any(pattern in host.lower() for pattern in rodent_patterns)

def write_host_distribution(f, sequences):
    """Write host distribution to summary file"""
    host_counts = {'human': 0, 'rodent': 0, 'other': 0, 'unknown': 0}
    other_hosts = set()
    
    for seq in sequences:
        record = seq['record']
        host_found = False
        
        for feature in record.features:
            if feature.type == "source" and 'host' in feature.qualifiers:
                host = feature.qualifiers['host'][0].lower()
                host_found = True
                
                if is_human_host(host):
                    host_counts['human'] += 1
                elif is_rodent_host(host):
                    host_counts['rodent'] += 1
                else:
                    host_counts['other'] += 1
                    other_hosts.add(host)
                break
        
        if not host_found:
            host_counts['unknown'] += 1
    
    f.write(f"Human hosts: {host_counts['human']}\n")
    f.write(f"Rodent hosts: {host_counts['rodent']}\n")
    f.write(f"Other hosts: {host_counts['other']}\n")
    f.write(f"Unknown hosts: {host_counts['unknown']}\n")
    
    if other_hosts:
        f.write("\nOther host types found:\n")
        for host in sorted(other_hosts):
            f.write(f"- {host}\n")
    f.write("\n")

def write_filtering_steps(f, sequences, genome_choice, completeness, host_choice, completeness_filtered, filtered_sequences):
    """Document the filtering steps and their impact"""
    f.write("\nFiltering Steps Summary\n")
    f.write("=====================\n\n")
    
    # Initial counts
    initial_counts = calculate_segment_counts(sequences)
    initial_locations = calculate_location_counts(sequences)
    
    f.write("1. Initial Dataset\n")
    f.write("----------------\n")
    f.write(f"Total sequences: {len(sequences)}\n")
    f.write(f"L segments: {initial_counts['L']}\n")
    f.write(f"S segments: {initial_counts['S']}\n")
    f.write(f"Unknown segments: {initial_counts['unknown']}\n\n")
    
    # Add detailed host information
    f.write("Initial Host Distribution\n")
    f.write("----------------------\n")
    host_counts = {}  # Dictionary to store counts for each unique host
    host_details = {'human': set(), 'rodent': set(), 'other': set()}
    no_host_count = 0
    
    for seq in sequences:
        record = seq['record']
        host_found = False
        
        for feature in record.features:
            if feature.type == "source" and 'host' in feature.qualifiers:
                host = feature.qualifiers['host'][0]
                host_found = True
                
                # Update count for this specific host
                if host not in host_counts:
                    host_counts[host] = 0
                host_counts[host] += 1
                
                if is_human_host(host):
                    host_details['human'].add(host.lower())
                elif is_rodent_host(host):
                    host_details['rodent'].add(host.lower())
                else:
                    host_details['other'].add(host.lower())
                break
        
        if not host_found:
            no_host_count += 1
    
    # Calculate totals for each category
    total_human = sum(host_counts[host] for host in host_counts if is_human_host(host))
    total_rodent = sum(host_counts[host] for host in host_counts if is_rodent_host(host))
    total_other = sum(host_counts[host] for host in host_counts 
                     if not is_human_host(host) and not is_rodent_host(host))
    
    # Write human hosts with counts
    f.write(f"\nHuman hosts found (Total: {total_human} sequences):\n")
    for host in sorted(host_details['human']):
        matching_hosts = [h for h in host_counts.keys() if h.lower() == host]
        for matching_host in matching_hosts:
            f.write(f"- {matching_host}: {host_counts[matching_host]} sequences\n")
    
    # Write rodent hosts with counts
    f.write(f"\nRodent hosts found (Total: {total_rodent} sequences):\n")
    for host in sorted(host_details['rodent']):
        matching_hosts = [h for h in host_counts.keys() if h.lower() == host]
        for matching_host in matching_hosts:
            f.write(f"- {matching_host}: {host_counts[matching_host]} sequences\n")
    
    # Write other hosts with counts
    if host_details['other']:
        f.write(f"\nOther hosts found (Total: {total_other} sequences):\n")
        for host in sorted(host_details['other']):
            matching_hosts = [h for h in host_counts.keys() if h.lower() == host]
            for matching_host in matching_hosts:
                f.write(f"- {matching_host}: {host_counts[matching_host]} sequences\n")
    
    f.write(f"\nSequences with no host information: {no_host_count}\n\n")
    
    write_geographical_distribution(f, initial_locations, "Initial Geographical Distribution")
    
    # Completeness filtering
    if genome_choice != '3':
        f.write("\n2. After Completeness Filtering\n")
        f.write("----------------------------\n")
        if genome_choice == '1':
            f.write("Applied filter: Complete genomes only (>99%)\n")
        else:
            f.write(f"Applied filter: Minimum completeness {completeness}%\n")
        
        completeness_counts = calculate_segment_counts(completeness_filtered)
        completeness_locations = calculate_location_counts(completeness_filtered)
        
        f.write(f"\nRemaining sequences: {len(completeness_filtered)}\n")
        f.write(f"L segments: {completeness_counts['L']}\n")
        f.write(f"S segments: {completeness_counts['S']}\n")
        f.write(f"Unknown segments: {completeness_counts['unknown']}\n\n")
        write_geographical_distribution(f, completeness_locations, "Geographical Distribution After Completeness Filter")
    
    # Host filtering
    if host_choice != '4':
        f.write("\n3. After Host Filtering\n")
        f.write("---------------------\n")
        host_type = {
            '1': 'Human only',
            '2': 'Rodent only',
            '3': 'Human and Rodent'
        }[host_choice]
        f.write(f"Applied filter: {host_type}\n")
        
        final_counts = calculate_segment_counts(filtered_sequences)
        final_locations = calculate_location_counts(filtered_sequences)
        
        f.write(f"\nFinal sequences: {len(filtered_sequences)}\n")
        f.write(f"L segments: {final_counts['L']}\n")
        f.write(f"S segments: {final_counts['S']}\n")
        f.write(f"Unknown segments: {final_counts['unknown']}\n\n")
        write_geographical_distribution(f, final_locations, "Final Geographical Distribution")

def filter_by_host(sequences, host_choice):
    """Filter sequences based on host"""
    filtered_sequences = []
    
    for seq in sequences:
        record = seq['record']
        include_sequence = False
        host_found = False
        
        for feature in record.features:
            if feature.type == "source" and 'host' in feature.qualifiers:
                host_found = True
                host = feature.qualifiers['host'][0].lower()
                
                if host_choice == '1':  # Human only
                    include_sequence = is_human_host(host)
                elif host_choice == '2':  # Rodent only
                    include_sequence = is_rodent_host(host)
                elif host_choice == '3':  # Both human and rodent
                    include_sequence = is_human_host(host) or is_rodent_host(host)
                elif host_choice == '4':  # No filter
                    include_sequence = True
                break
        
        if host_choice == '4' or include_sequence:
            filtered_sequences.append(seq)
    
    return filtered_sequences

def filter_by_metadata(sequences, metadata_choice):
    """Filter sequences based on metadata completeness"""
    filtered_sequences = []
    
    for seq in sequences:
        record = seq['record']
        metadata = get_metadata(record)  # This gives us accession_location_date
        
        include_sequence = True
        if metadata_choice == '1':  # Known location only
            if 'UnknownLoc' in metadata:
                include_sequence = False
        elif metadata_choice == '2':  # Known date only
            if 'UnknownDate' in metadata:
                include_sequence = False
        elif metadata_choice == '3':  # Both known
            if 'UnknownLoc' in metadata or 'UnknownDate' in metadata:
                include_sequence = False
        # metadata_choice == '4' means no filter, include all sequences
        
        if include_sequence:
            filtered_sequences.append(seq)
    
    return filtered_sequences

def filter_by_country(sequences, countries):
    """Filter sequences based on specified countries"""
    if not countries:  # If no countries specified, return all sequences
        return sequences
        
    filtered_sequences = []
    # Convert input countries to standardized format
    standardized_countries = [get_standardized_country_name(c) for c in countries]
    
    for seq in sequences:
        record = seq['record']
        include_sequence = False
        
        for feature in record.features:
            if feature.type == "source" and 'geo_loc_name' in feature.qualifiers:
                location = feature.qualifiers['geo_loc_name'][0]
                if 'missing' not in location.lower():
                    country = location.split(':')[0].strip()
                    standardized_country = get_standardized_country_name(country)
                    if standardized_country in standardized_countries:
                        include_sequence = True
                        break
        
        if include_sequence:
            filtered_sequences.append(seq)
    
    return filtered_sequences

def write_metadata_filtering_summary(f, sequences, metadata_filtered_sequences, metadata_choice, written_counts):
    """Write metadata filtering summary to the report"""
    # Only write metadata filtering summary if actual filtering was applied
    if metadata_choice != '4':
        f.write("\n4. After Metadata Filtering\n")
        f.write("------------------------\n")
        
        # Add description of applied filter
        filter_descriptions = {
            '1': 'Applied filter: Known location only',
            '2': 'Applied filter: Known date only',
            '3': 'Applied filter: Both known location and date'
        }
        f.write(f"{filter_descriptions[metadata_choice]}\n")
        
        # Count remaining sequences after filtering
        final_counts = calculate_segment_counts(metadata_filtered_sequences)
        final_locations = calculate_location_counts(metadata_filtered_sequences)
        
        f.write(f"\nAfter metadata filtering:\n")
        f.write(f"Final sequences: {len(metadata_filtered_sequences)}\n")
        f.write(f"L segments: {final_counts['L']}\n")
        f.write(f"S segments: {final_counts['S']}\n")
        f.write(f"Unknown segments: {final_counts['unknown']}\n\n")
        
        write_geographical_distribution(f, final_locations, "Final Geographical Distribution After Metadata Filter")

def write_country_filtering_summary(f, sequences, country_filtered_sequences, countries, written_counts):
    """Write country filtering summary and final output files section"""
    if countries:  # Only write if country filtering was applied
        f.write("\n5. After Country Filtering\n")
        f.write("----------------------\n")
        f.write(f"Applied filter: Selected countries: {countries}\n")
        
        final_counts = calculate_segment_counts(country_filtered_sequences)
        final_locations = calculate_location_counts(country_filtered_sequences)
        
        f.write(f"\nAfter country filtering:\n")
        f.write(f"Final sequences: {len(country_filtered_sequences)}\n")
        f.write(f"L segments: {final_counts['L']}\n")
        f.write(f"S segments: {final_counts['S']}\n")
        f.write(f"Unknown segments: {final_counts['unknown']}\n\n")
        
        write_geographical_distribution(f, final_locations, "Final Geographical Distribution After Country Filter")
    
    # Always write the Output Files section at the very end
    f.write("\nOutput Files\n")
    f.write("------------\n")
    f.write(f"L segments: FASTA/L_segment/lassa_l_segments.fasta ({written_counts['L']} sequences)\n")
    f.write(f"S segments: FASTA/S_segment/lassa_s_segments.fasta ({written_counts['S']} sequences)\n")
    if written_counts['unknown'] > 0:
        f.write(f"Unknown segments: FASTA/unknown_segment/lassa_unknown_segments.fasta ({written_counts['unknown']} sequences)\n")

def create_phylogeny_directories(output_dir):
    """Create all necessary directories for phylogenetic analysis"""
    phylogeny_dir = os.path.join(output_dir, "Phylogeny")
    
    # Create main phylogeny subdirectories (removed Recombination)
    subdirs = ["FASTA", "MSA", "TrimAl", "Tree"]
    for subdir in subdirs:
        # Create L and S segment subdirectories in each main directory
        l_dir = os.path.join(phylogeny_dir, subdir, "L_segment")
        s_dir = os.path.join(phylogeny_dir, subdir, "S_segment")
        os.makedirs(l_dir, exist_ok=True)
        os.makedirs(s_dir, exist_ok=True)
    
    return phylogeny_dir

def copy_consensus_sequence(consensus_file, segment, output_dir):
    """Copy consensus sequence(s) to the appropriate directory and format their headers.
    
    Args:
        consensus_file (str): Path to consensus sequence file in FASTA format
        segment (str): Segment type ('L' or 'S')
        output_dir (str): Output directory path
    """
    from Bio import SeqIO
    import shutil
    
    if not consensus_file:
        return
        
    # Create consensus directory if it doesn't exist
    segment_dir = os.path.join(output_dir, "FASTA", f"{segment}_segment")
    os.makedirs(segment_dir, exist_ok=True)
    
    # Read all sequences from the consensus file
    records = []
    for record in SeqIO.parse(consensus_file, "fasta"):
        # Format the header: OriginalID_L_segment or OriginalID_S_segment
        record.id = f"{record.id}_{segment}_segment"
        record.description = ""
        records.append(record)
    
    # Write all sequences to consensus.fasta in the segment directory
    output_file = os.path.join(segment_dir, "consensus.fasta")
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")
    
    print(f"Added {len(records)} consensus sequence(s) for {segment} segment")

def concatenate_fasta_files(input_dir, phylogeny_dir, segment):
    """Concatenate all FASTA files in a segment directory and remove duplicates"""
    from Bio import SeqIO
    from collections import defaultdict
    
    # Initialize dictionaries to store sequences and track duplicates
    sequences = {}
    duplicate_count = 0
    duplicate_headers = defaultdict(list)
    
    # Define reference lengths for L and S segments
    reference_lengths = {
        'L': 7200,  # ~7.2kb for L segment
        'S': 3400   # ~3.4kb for S segment
    }
    min_length = reference_lengths[segment] * 0.10  # 10% of reference length
    
    # Read from the original FASTA directory
    segment_dir = os.path.join(input_dir, "FASTA", f"{segment}_segment")
    
    # Read all FASTA files in the directory
    for fasta_file in ["lassa_" + segment.lower() + "_segments.fasta", 
                      "reference.fasta", "outgroup.fasta", "consensus.fasta"]:  # Added consensus.fasta
        file_path = os.path.join(segment_dir, fasta_file)
        if os.path.exists(file_path):
            for record in SeqIO.parse(file_path, "fasta"):
                # Skip sequences shorter than 10% of reference length for phylogeny
                if len(record.seq) < min_length:
                    print(f"Skipping {record.id} (length: {len(record.seq)}bp) - shorter than minimum required length ({min_length}bp) for phylogenetic analysis")
                    continue
                    
                # Check if this header already exists
                if record.id in sequences:
                    duplicate_count += 1
                    duplicate_headers[record.id].append(fasta_file)
                else:
                    sequences[record.id] = record
    
    # Write only to the segment subdirectory in FASTA
    segment_output_dir = os.path.join(phylogeny_dir, "FASTA", f"{segment}_segment")
    os.makedirs(segment_output_dir, exist_ok=True)
    
    # Write concatenated sequences to new file
    output_file = os.path.join(segment_output_dir, f"all_{segment.lower()}_segments.fasta")
    with open(output_file, "w") as f:
        SeqIO.write(sequences.values(), f, "fasta")
    
    # Create metadata file in the same directory
    metadata_output = os.path.join(segment_output_dir, f"{segment.lower()}_metadata.txt")
    create_figtree_metadata(output_file, metadata_output)
    print(f"Created FigTree metadata file: {metadata_output}")
    
    return len(sequences)  # Return number of unique sequences

def download_and_write_special_sequences(output_dir):
    """Download reference and outgroup sequences and write to appropriate files"""
    from Bio import Entrez
    from Bio import SeqIO
    from io import StringIO
    Entrez.email = "your.email@example.com"
    
    for segment in ['L', 'S']:
        segment_dir = os.path.join(output_dir, "FASTA", f"{segment}_segment")
        
        # Download and format reference sequence
        ref_handle = Entrez.efetch(db="nucleotide", 
                                 id=REFERENCE_SEQUENCES[segment]['id'], 
                                 rettype="fasta", 
                                 retmode="text")
        
        # Parse and format reference sequence header
        ref_record = SeqIO.read(StringIO(ref_handle.read()), "fasta")
        ref_record.id = f"{REFERENCE_SEQUENCES[segment]['id']}_SierraLeone_Human_Reference_Unknown"  # Moved date to end
        ref_record.description = ""
        
        # Write reference sequence with formatted header
        with open(os.path.join(segment_dir, "reference.fasta"), "w") as f:
            SeqIO.write(ref_record, f, "fasta")
        
        # Download and format outgroup sequence
        out_handle = Entrez.efetch(db="nucleotide", 
                                 id=PINNEO_SEQUENCES[segment]['id'], 
                                 rettype="fasta", 
                                 retmode="text")
        
        # Parse and format outgroup sequence header
        out_record = SeqIO.read(StringIO(out_handle.read()), "fasta")
        out_record.id = f"{PINNEO_SEQUENCES[segment]['id']}_Nigeria_Lassa_Human_Pinneo_outgroup_{PINNEO_SEQUENCES[segment]['date']}"  # Changed SierraLeone to Nigeria
        out_record.description = ""
        
        # Write outgroup sequence with formatted header
        with open(os.path.join(segment_dir, "outgroup.fasta"), "w") as f:
            SeqIO.write(out_record, f, "fasta")

def get_lineages_dir():
    """Get the path to the lineages directory from the LassaSeq package installation"""
    import os
    import sys
    import site
    
    def find_lassaseq_root():
        # First try: Look in the conda environment's site-packages
        conda_prefix = os.environ.get('CONDA_PREFIX')
        if conda_prefix:
            python_version = f"python{sys.version_info.major}.{sys.version_info.minor}"
            site_packages = os.path.join(conda_prefix, 'lib', python_version, 'site-packages')
            lassaseq_dir = os.path.join(site_packages, 'lassaseq')
            if os.path.exists(os.path.join(lassaseq_dir, 'lineages')):
                return lassaseq_dir
        
        # Second try: Look in site-packages
        for site_dir in site.getsitepackages():
            lassaseq_dir = os.path.join(site_dir, 'lassaseq')
            if os.path.exists(os.path.join(lassaseq_dir, 'lineages')):
                return lassaseq_dir
        
        # Third try: Check if we're in development mode
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if os.path.exists(os.path.join(script_dir, 'lineages')):
            return script_dir
        
        # Fourth try: Look in current directory
        cwd = os.getcwd()
        if os.path.exists(os.path.join(cwd, 'lineages')):
            return cwd
        
        # Fifth try: Look in LassaSeq directory if we're in it
        if os.path.basename(cwd) == 'LassaSeq':
            lassaseq_dir = os.path.join(cwd, 'lassaseq')
            if os.path.exists(os.path.join(lassaseq_dir, 'lineages')):
                return lassaseq_dir
        
        return None
    
    # Find LassaSeq root directory
    lassaseq_root = find_lassaseq_root()
    if lassaseq_root is None:
        print(f"Error: Could not find LassaSeq installation directory")
        print("Please ensure you are running from within the LassaSeq directory or that LassaSeq is properly installed")
        sys.exit(1)
    
    # Get the lineages directory path
    lineages_dir = os.path.join(lassaseq_root, 'lineages')
    
    # Verify the lineages directory exists
    if not os.path.exists(lineages_dir):
        print(f"Error: Lineages directory not found at {lineages_dir}")
        print("Please ensure the lineages directory is present in your LassaSeq installation")
        sys.exit(1)
    
    return lineages_dir

def create_figtree_metadata(trimmed_fasta, output_file, segment=None):
    """Create metadata file for FigTree visualization from trimmed alignment headers."""
    from Bio import SeqIO
    import os
    
    # Keep track of processed sequences to avoid duplicates
    processed_sequences = set()
    
    # Load lineage information from files
    lineages = {}
    sublineages = {}
    subsublineages = {}
    
    # Get the lineages directory path
    lineages_dir = get_lineages_dir()
    
    # Determine if this is L or S segment based on explicit arg, protein folder, or path
    if segment is None:
        path_u = output_file.replace('\\', '/')
        if 'L_segment' in path_u or '/L/' in path_u or path_u.rstrip('/').endswith('/L'):
            segment = 'l'
        elif 'S_segment' in path_u or '/S/' in path_u or path_u.rstrip('/').endswith('/S'):
            segment = 's'
        else:
            # Protein folders: L/Z → L segment; NP/GPC → S segment
            for prot, info in LASSA_PROTEIN_REGIONS.items():
                if f'/{prot}/' in path_u or f'/{prot.lower()}/' in path_u.lower():
                    segment = info['segment'].lower()
                    break
            else:
                segment = 'l' if '_l_' in path_u.lower() else 's'
    else:
        segment = segment.lower()
    
    # Get the path to the lineages file
    lineage_file = os.path.join(lineages_dir, f'{segment}_lineages.txt')
    
    # Read lineage information if file exists
    if os.path.exists(lineage_file):
        with open(lineage_file) as f:
            header = [h.strip() for h in next(f).strip().split('\t')]
            header_l = [h.lower() for h in header]
            acc_idx = next((i for i, col in enumerate(header_l) if 'acc' in col), 0)
            # Prefer exact-ish column names so "Sublineage" is not mistaken for Lineages
            lin_idx = next((i for i, col in enumerate(header_l)
                            if col.startswith('lineage') or col == 'lineages'), 1)
            sub_idx = next((i for i, col in enumerate(header_l)
                            if col == 'sublineage' or (col.startswith('sub') and 'subsub' not in col)), 2)
            subsub_idx = next((i for i, col in enumerate(header_l) if 'subsub' in col), None)
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    accession = parts[acc_idx].strip()
                    lineages[accession] = parts[lin_idx].strip()
                    sublineages[accession] = parts[sub_idx].strip()
                    if subsub_idx is not None and len(parts) > subsub_idx:
                        subsublineages[accession] = parts[subsub_idx].strip() or "None"
                    else:
                        subsublineages[accession] = "None"
    
    # Write header line first
    with open(output_file, 'w') as f:
        f.write("Taxon\tLocation\tLocation2\tHost\tDate\tLineage\tSublineage\tSubsublineage\n")
        
        # Parse FASTA file and extract metadata from headers
        for record in SeqIO.parse(trimmed_fasta, "fasta"):
            # Get the full header as taxon
            taxon = record.id
            
            # Skip if we've already processed this sequence
            if record.id in processed_sequences:
                continue
                
            # Skip consensus sequences (identified by _segment suffix)
            if record.id.endswith('_segment'):
                continue
            
            # Split the header into parts
            parts = record.id.split('_')
            
            processed_sequences.add(record.id)
            
            # Get accession from the first part and clean it
            accession = parts[0]
            # Handle special cases
            if "Reference" in parts:
                accession = "NC_004297.1" if segment == 'l' else "NC_004296.1"
            elif "outgroup" in parts:
                accession = "KM822127.1" if segment == 'l' else "KM822128.1"
            
            # Get lineage information
            lineage = lineages.get(accession, "Unknown")
            sublineage = sublineages.get(accession, "None")
            subsublineage = subsublineages.get(accession, "None")
            
            # Handle reference sequences specially
            if "Reference" in parts:
                f.write(f"{taxon}\tSierraLeone\tUnknown\tHuman\tUnknown\t{lineage}\t{sublineage}\t{subsublineage}\n")
            elif "outgroup" in parts:
                f.write(f"{taxon}\tNigeria\tLassa\tHuman\t1969.000\t{lineage}\t{sublineage}\t{subsublineage}\n")
            else:
                location = parts[1] if len(parts) > 1 else "Unknown"
                location2 = parts[2] if len(parts) > 2 else "Unknown"
                host = parts[3] if len(parts) > 3 else "Unknown"
                date = parts[4] if len(parts) > 4 else "Unknown"
                f.write(f"{taxon}\t{location}\t{location2}\t{host}\t{date}\t{lineage}\t{sublineage}\t{subsublineage}\n")

def perform_phylogenetic_analysis(phylogeny_dir, segment, do_tree=True):
    """Perform alignment, trimming and optionally tree building for a segment"""
    from Bio import SeqIO
    
    # Define reference lengths and minimum required length
    reference_lengths = {
        'L': 7200,  # ~7.2kb for L segment
        'S': 3400   # ~3.4kb for S segment
    }
    min_length = reference_lengths[segment] * 0.10  # 10% of reference length
    
    # Setup directories
    fasta_dir = os.path.join(phylogeny_dir, "FASTA", f"{segment}_segment")
    msa_dir = os.path.join(phylogeny_dir, "MSA", f"{segment}_segment")
    trimal_dir = os.path.join(phylogeny_dir, "TrimAl", f"{segment}_segment")
    tree_dir = os.path.join(phylogeny_dir, "Tree", f"{segment}_segment")
    
    input_fasta = os.path.join(fasta_dir, f"all_{segment.lower()}_segments.fasta")
    
    # First filter sequences by length
    filtered_records = []
    skipped_count = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        if len(record.seq) >= min_length:
            filtered_records.append(record)
        else:
            skipped_count += 1
            print(f"Skipping {record.id} (length: {len(record.seq)}bp) - shorter than minimum required length ({min_length}bp) for phylogenetic analysis")
    
    if skipped_count > 0:
        print(f"\nSkipped {skipped_count} sequences that were shorter than {min_length}bp (10% of reference length)")
    
    if not filtered_records:
        print(f"Error: No sequences remaining after length filtering for {segment} segment")
        return
    
    # Write filtered sequences to MSA directory
    msa_input = os.path.join(msa_dir, f"{segment.lower()}_sequences.fasta")
    with open(msa_input, "w") as f:
        SeqIO.write(filtered_records, f, "fasta")
    
    aligned_output = os.path.join(msa_dir, f"{segment.lower()}_aligned.fasta")
    temp_aligned = os.path.join(msa_dir, "temp_aligned.fasta")
    
    # Store original headers
    original_headers = {}
    for record in filtered_records:
        original_headers[record.id] = record.id
    
    # Perform MAFFT alignment
    print(f"\nAligning {segment} segment sequences using MAFFT...")
    
    mafft_cmd = [
        "mafft",
        "--localpair",
        "--maxiterate", "1000",
        "--ep", "0.123",
        "--op", "1.53",
        "--reorder",
        "--adjustdirection",
        "--thread", "-1",
        msa_input
    ]
    
    try:
        with open(temp_aligned, 'w') as outfile:
            process = subprocess.Popen(
                mafft_cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True
            )
            _, stderr = process.communicate()
            
            if process.returncode != 0:
                print(f"MAFFT Error: {stderr}")
                raise subprocess.CalledProcessError(process.returncode, mafft_cmd)
        
        # Restore original headers while keeping aligned sequences
        with open(aligned_output, 'w') as outfile:
            for record in SeqIO.parse(temp_aligned, "fasta"):
                # Remove _R_ prefix if present and get original header
                orig_id = original_headers[record.id.replace('_R_', '')]
                outfile.write(f">{orig_id}\n{str(record.seq)}\n")
        
        # Remove temporary files
        os.remove(temp_aligned)
        os.remove(msa_input)
        
        print(f"MAFFT alignment completed for {segment} segment")
        
        # Copy aligned file to TrimAl directory and perform trimming
        trimal_input = os.path.join(trimal_dir, f"{segment.lower()}_aligned.fasta")
        trimal_output = os.path.join(trimal_dir, f"{segment.lower()}_trimmed.fasta")
        shutil.copy2(aligned_output, trimal_input)
        
        print(f"\nTrimming {segment} segment alignment using TrimAl...")
        
        trimal_cmd = [
            "trimal",
            "-in", trimal_input,
            "-out", trimal_output,
            "-automated1"
        ]
        
        process = subprocess.run(trimal_cmd, capture_output=True, text=True)
        if process.returncode != 0:
            print(f"TrimAl Error: {process.stderr}")
            raise subprocess.CalledProcessError(process.returncode, trimal_cmd)
            
        print(f"TrimAl trimming completed for {segment} segment")
        
        # Remove input alignment file, keeping only trimmed output
        os.remove(trimal_input)

        if not do_tree:
            return
        
        # Copy trimmed alignment to Tree directory
        tree_input = os.path.join(tree_dir, f"{segment.lower()}_trimmed.fasta")
        shutil.copy2(trimal_output, tree_input)
        
        # Copy metadata file from FASTA directory to Tree directory
        fasta_metadata = os.path.join(fasta_dir, f"{segment.lower()}_metadata.txt")
        tree_metadata = os.path.join(tree_dir, f"{segment.lower()}_metadata.txt")
        shutil.copy2(fasta_metadata, tree_metadata)
        print(f"Copied metadata file to Tree folder: {tree_metadata}")
        
        iqtree_bin = _resolve_iqtree_cmd() or "iqtree2"
        print(f"\nBuilding phylogenetic tree for {segment} segment using {iqtree_bin}...")
        
        # Change to tree directory before running IQ-TREE
        current_dir = os.getcwd()
        os.chdir(tree_dir)
        
        try:
            # Run IQ-TREE with just the filename since we're in the correct directory
            iqtree_cmd = [
                iqtree_bin,
                "-s", f"{segment.lower()}_trimmed.fasta",
                "-nt", "AUTO",
                "-bb", "10000",
                "-m", "TEST"
            ]
            
            process = subprocess.run(iqtree_cmd, capture_output=True, text=True)
            if process.returncode != 0:
                print(f"IQ-TREE Error: {process.stderr}")
                raise subprocess.CalledProcessError(process.returncode, iqtree_cmd)
                
            print(f"IQ-TREE analysis completed for {segment} segment")
            
            # Clean up intermediate files, keeping only essential outputs
            if os.path.exists(tree_dir):
                for file in os.listdir('.'):  # Use '.' since we're already in tree_dir
                    if file.endswith(('.bionj', '.ckp.gz', '.mldist', '.log')):
                        try:
                            os.remove(file)
                        except OSError as e:
                            print(f"Warning: Could not remove {file}: {e}")
                    
        finally:
            # Always return to the original directory
            os.chdir(current_dir)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {str(e)}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        raise


# ---------------------------------------------------------------------------
# Protein (CDS) alignment pipeline (minimap2 → translate → MAFFT → pal2nal)
# ---------------------------------------------------------------------------

def normalize_protein_list(proteins_str):
    """Parse comma-separated protein names/numbers into canonical names."""
    if not proteins_str:
        return []
    raw = [x.strip().upper() for x in proteins_str.split(',') if x.strip()]
    return [PROTEIN_NUM.get(p, p) for p in raw]


def check_protein_cds_dependencies():
    """Require minimap2, mafft, trimal, and pal2nal.pl for protein alignment."""
    for tool in ['minimap2', 'mafft', 'trimal']:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            print(f"Error: {tool} not found. Required for protein alignment.")
            sys.exit(1)
    if subprocess.run(['which', 'pal2nal.pl'], capture_output=True).returncode != 0:
        print("Error: pal2nal.pl not found. Install with e.g. conda install -c bioconda pal2nal")
        sys.exit(1)


def _resolve_iqtree_cmd():
    """Prefer iqtree2, then iqtree / iqtree3."""
    for cmd in ('iqtree2', 'iqtree', 'iqtree3'):
        if subprocess.run(['which', cmd], capture_output=True).returncode == 0:
            return cmd
    return None


def _fetch_refseq_fasta(accession):
    """Fetch RefSeq nucleotide FASTA from NCBI."""
    try:
        from urllib.request import urlopen
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        )
        with urlopen(url, timeout=30) as r:
            return r.read().decode()
    except Exception as e:
        raise RuntimeError(f"Failed to fetch RefSeq {accession}: {e}")


def _parse_cigar(cigar_str):
    for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_str):
        yield m.group(2), int(m.group(1))


def _extract_region_cds(seq, ref_start_1, ref_end_1, align_ref_start_1, cigar_str):
    """Extract query bases that align to a 1-based inclusive RefSeq interval."""
    out, ref_pos, query_pos = [], align_ref_start_1, 0
    for op, length in _parse_cigar(cigar_str):
        if op in "M=X":
            for _ in range(length):
                if ref_start_1 <= ref_pos <= ref_end_1:
                    out.append(seq[query_pos])
                ref_pos += 1
                query_pos += 1
        elif op == "I":
            if ref_start_1 <= ref_pos - 1 <= ref_end_1:
                for i in range(length):
                    out.append(seq[query_pos + i])
            query_pos += length
        elif op in "D N":
            ref_pos += length
        elif op in "S H P":
            query_pos += length
    return "".join(out)


def _revcomp(seq):
    complement = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
    return seq.translate(complement)[::-1]


def _extract_cds_from_sam(sam_path, ref_start, ref_end, query_ids, out_path,
                          min_cds_fraction=0.5, strand='+'):
    """Extract CDS nucleotides from a minimap2 SAM; reverse-complement if strand='-'."""
    ref_len = ref_end - ref_start + 1
    min_len = int(ref_len * min_cds_fraction)
    query_set = set(query_ids) if query_ids else None
    seqs = OrderedDict()
    best_extracted_len = {}
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            qname, flag, rname, pos_s, _, cigar, _, _, _, seq = parts[:10]
            if int(flag) & 0x900 or not seq or seq == "*":
                continue
            pos = int(pos_s)
            if pos > ref_end:
                continue
            if query_set is not None and qname not in query_set:
                continue
            extracted = _extract_region_cds(seq, ref_start, ref_end, pos, cigar)
            if strand == '-':
                extracted = _revcomp(extracted)
            n_ext = len(extracted)
            best_extracted_len[qname] = max(best_extracted_len.get(qname, 0), n_ext)
            if n_ext < min_len:
                continue
            if qname not in seqs or n_ext > len(seqs[qname]):
                seqs[qname] = extracted
    with open(out_path, "w") as f:
        for qid in (query_ids or []):
            if qid in seqs:
                s = seqs[qid]
                f.write(f">{qid}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + "\n")
    kept = set(seqs.keys())
    dropped_ids = [q for q in (query_ids or []) if q not in kept]
    dropped_best_bp = {q: best_extracted_len.get(q, 0) for q in dropped_ids}
    return len(seqs), len(query_ids or []), dropped_ids, ref_len, dropped_best_bp


def _translate_cds_to_protein(rec):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    s = rec.seq
    if len(s) % 3 != 0:
        s = s[:(len(s) // 3) * 3]
    aa = str(s.translate(table=1, to_stop=False, cds=False))
    return SeqRecord(Seq(aa), id=rec.id, description="")


def _translate_cds_fasta(cds_path, protein_path):
    records = list(SeqIO.parse(cds_path, "fasta"))
    if not records:
        open(protein_path, "w").close()
        return
    proteins = [_translate_cds_to_protein(rec) for rec in records]
    with open(protein_path, "w") as f:
        SeqIO.write(proteins, f, "fasta")


def _drop_proteins_with_internal_stops(cds_path, protein_path, protein_label):
    """Remove CDS/protein pairs with internal stop codons before MAFFT."""
    proteins = list(SeqIO.parse(protein_path, "fasta"))
    if not proteins:
        return 0
    bad_ids = []
    for rec in proteins:
        aa = str(rec.seq).upper().replace("-", "")
        body = aa[:-1] if aa.endswith("*") else aa
        if "*" in body:
            bad_ids.append(rec.id)
    if not bad_ids:
        return 0
    bad_set = set(bad_ids)
    for eid in bad_ids:
        shown = eid if len(eid) <= 80 else eid[:80] + "..."
        print(f"  {protein_label}: excluding before alignment (internal stop codon): {shown}")
    kept_prot = [r for r in proteins if r.id not in bad_set]
    kept_cds = [r for r in SeqIO.parse(cds_path, "fasta") if r.id not in bad_set]
    with open(protein_path, "w") as f:
        SeqIO.write(kept_prot, f, "fasta")
    with open(cds_path, "w") as f:
        SeqIO.write(kept_cds, f, "fasta")
    print(
        f"  {protein_label}: removed {len(bad_ids)} sequence(s) with internal stop(s) "
        f"before MAFFT ({len(kept_prot)} remain)."
    )
    return len(bad_ids)


def _parse_pal2nal_failed_id(err):
    if not err:
        return None
    for line in err.splitlines():
        s = line.strip()
        if s.startswith(">"):
            return s[1:].strip()
    m = re.search(r">\s*([^\s\n\r>]+)", err)
    return m.group(1).strip() if m else None


def _backtranslate_pal2nal(protein_aln_path, cds_path, out_path):
    """Back-translate protein alignment to codon alignment using PAL2NAL."""
    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    cds_by_id = {}
    for rec in SeqIO.parse(cds_path, "fasta"):
        cds_by_id[rec.id] = rec.seq
        cds_by_id[rec.id.replace("/", "_")] = rec.seq
        cds_by_id[rec.id.replace("_", "/")] = rec.seq
    aln = AlignIO.read(protein_aln_path, "fasta")
    cds_ordered = []
    for rec in aln:
        seq = (
            cds_by_id.get(rec.id)
            or cds_by_id.get(rec.id.replace("/", "_"))
            or cds_by_id.get(rec.id.replace("_", "/"))
        )
        if seq is None:
            return False, f"CDS missing for id: {rec.id}", None
        cds_ordered.append(SeqRecord(seq, id=rec.id, description=""))
    fd, tmp_cds = tempfile.mkstemp(suffix=".fasta")
    os.close(fd)
    try:
        with open(tmp_cds, "w") as f:
            SeqIO.write(cds_ordered, f, "fasta")
        ret = subprocess.run(
            ["pal2nal.pl", protein_aln_path, tmp_cds, "-output", "fasta"],
            capture_output=True, text=True, timeout=600,
        )
        err_text = ((ret.stderr or "").strip() + "\n" + (ret.stdout or "").strip()).strip()
        if ret.returncode != 0:
            err = err_text or f"exit code {ret.returncode}"
            return False, err, _parse_pal2nal_failed_id(err)
        with open(out_path, "w") as f:
            f.write(ret.stdout)
        n_out = len(list(SeqIO.parse(out_path, "fasta")))
        if n_out == 0:
            err = err_text or "PAL2NAL produced 0 sequences (protein/CDS mismatch?)"
            return False, err, _parse_pal2nal_failed_id(err)
        return True, None, None
    finally:
        try:
            os.unlink(tmp_cds)
        except OSError:
            pass


def _collect_segment_fasta_for_protein(output_dir, segment):
    """Build a deduplicated FASTA for one segment (download + ref + outgroup + consensus)."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    reference_lengths = {'L': 7200, 'S': 3400}
    min_length = reference_lengths[segment] * 0.10
    sequences = OrderedDict()
    segment_dir = os.path.join(output_dir, "FASTA", f"{segment}_segment")
    for fasta_file in [
        f"lassa_{segment.lower()}_segments.fasta",
        "reference.fasta",
        "outgroup.fasta",
        "consensus.fasta",
    ]:
        file_path = os.path.join(segment_dir, fasta_file)
        if not os.path.exists(file_path):
            continue
        for record in SeqIO.parse(file_path, "fasta"):
            if len(record.seq) < min_length:
                continue
            if record.id not in sequences:
                sequences[record.id] = SeqRecord(
                    Seq(str(record.seq)), id=record.id, description=""
                )
    return list(sequences.values())


def run_protein_cds_pipeline(output_dir, proteins, min_cds_fraction=0.5, threads=1):
    """
    Extract Lassa protein CDS from L/S segment genomes, align at protein level,
    back-translate with pal2nal, and trim. Writes Alignment/MAFFT/<protein>/ and
    Alignment/Trimmed/<protein>/.
    """
    check_protein_cds_dependencies()
    alignment_dir = os.path.join(output_dir, "Alignment")
    work_dir = os.path.join(alignment_dir, "_work")
    os.makedirs(work_dir, exist_ok=True)

    # Prepare per-segment input FASTAs once
    segment_inputs = {}
    for segment in ('L', 'S'):
        records = _collect_segment_fasta_for_protein(output_dir, segment)
        if not records:
            print(f"Warning: no sequences found for {segment} segment; "
                  f"proteins on this segment will be skipped.")
            continue
        seg_fa = os.path.join(work_dir, f"{segment}_input.fasta")
        with open(seg_fa, "w") as f:
            SeqIO.write(records, f, "fasta")
        segment_inputs[segment] = (seg_fa, [r.id for r in records])
        print(f"Protein pipeline input ({segment} segment): {len(records)} sequences")

    # Map each segment once to its RefSeq
    segment_sam = {}
    for segment, (seg_fa, query_ids) in segment_inputs.items():
        ref_id = REFERENCE_SEQUENCES[segment]['id']
        ref_fa = os.path.join(work_dir, f"{segment}_ref.fa")
        try:
            with open(ref_fa, "w") as f:
                f.write(_fetch_refseq_fasta(ref_id))
        except RuntimeError as e:
            print(f"Warning: could not fetch {ref_id} ({e}); using local reference.fasta if present")
            local_ref = os.path.join(output_dir, "FASTA", f"{segment}_segment", "reference.fasta")
            if not os.path.exists(local_ref):
                print(f"Error: no reference available for {segment} segment")
                continue
            shutil.copy2(local_ref, ref_fa)
        sam_path = os.path.join(work_dir, f"{segment}_aln.sam")
        print(f"Mapping {segment} segment sequences to {ref_id} with minimap2...")
        with open(sam_path, "w") as f:
            subprocess.run(
                ["minimap2", "-a", "-t", str(threads), ref_fa, seg_fa],
                stdout=f, stderr=subprocess.DEVNULL,
            )
        segment_sam[segment] = (sam_path, query_ids)

    extraction_counts = {}
    for protein in proteins:
        info = LASSA_PROTEIN_REGIONS.get(protein)
        if info is None:
            print(f"Warning: unknown protein '{protein}', skipping")
            continue
        segment = info['segment']
        if segment not in segment_sam:
            print(f"  {protein}: skipped (no {segment} segment input)")
            continue
        sam_path, query_ids = segment_sam[segment]
        protdir = os.path.join(work_dir, protein)
        os.makedirs(protdir, exist_ok=True)
        cds_path = os.path.join(protdir, "cds.fasta")
        n_passed, n_input, dropped_ids, ref_len, dropped_best_bp = _extract_cds_from_sam(
            sam_path, info['start'], info['end'], query_ids, cds_path,
            min_cds_fraction=min_cds_fraction, strand=info['strand'],
        )
        extraction_counts[protein] = (n_passed, n_input, dropped_ids, ref_len, dropped_best_bp)
        print(
            f"  {protein}: extracted CDS for {n_passed}/{n_input} sequences "
            f"(dropped {n_input - n_passed} with coverage < {min_cds_fraction * 100:.0f}% of RefSeq CDS)"
        )
        if dropped_ids:
            for did in dropped_ids[:10]:
                best = dropped_best_bp.get(did, 0)
                pct = (100.0 * best / ref_len) if ref_len else 0
                print(f"    dropped: {did} → {best} bp ({pct:.1f}% of ref CDS)")
            if len(dropped_ids) > 10:
                print(f"    ... and {len(dropped_ids) - 10} more")

    # Align each protein
    for protein in proteins:
        protdir = os.path.join(work_dir, protein)
        combined_cds = os.path.join(protdir, "cds.fasta")
        if not os.path.isfile(combined_cds) or os.path.getsize(combined_cds) == 0:
            continue
        protein_fasta = os.path.join(protdir, "protein.fasta")
        _translate_cds_fasta(combined_cds, protein_fasta)
        _drop_proteins_with_internal_stops(combined_cds, protein_fasta, protein)
        if not list(SeqIO.parse(protein_fasta, "fasta")):
            print(f"  {protein}: no sequences left after stop-codon filter; skipping.")
            continue

        protein_aln = os.path.join(protdir, "protein_aln.fasta")
        print(f"  Aligning {protein} proteins with MAFFT...")
        with open(protein_aln, "w") as f:
            r = subprocess.run(
                ["mafft", "--auto", "--thread", str(min(threads, 32)), protein_fasta],
                stdout=f, stderr=subprocess.PIPE, text=True,
            )
            if r.returncode != 0:
                print(f"MAFFT failed for {protein}: {r.stderr}")
                continue

        cds_aligned = os.path.join(protdir, "cds_aligned.fasta")
        exclude_for_pal2nal = set()
        pal2nal_ok = False
        for _ in range(50):
            cds_input = combined_cds
            aln_input = protein_aln
            if exclude_for_pal2nal:
                from Bio import AlignIO
                from Bio.Align import MultipleSeqAlignment

                def _excluded(eid):
                    return eid in exclude_for_pal2nal

                aln = AlignIO.read(protein_aln, "fasta")
                aln_filtered = [r for r in aln if not _excluded(r.id)]
                if not aln_filtered:
                    break
                cds_by_id = {rec.id: rec for rec in SeqIO.parse(combined_cds, "fasta")}
                cds_records = [cds_by_id[r.id] for r in aln_filtered if r.id in cds_by_id]
                if len(cds_records) != len(aln_filtered):
                    print(f"PAL2NAL retry: ID mismatch for {protein}, skipping retry.")
                    break
                fd_a = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
                fd_c = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
                with open(fd_a.name, "w") as f:
                    AlignIO.write(MultipleSeqAlignment(aln_filtered), f, "fasta")
                with open(fd_c.name, "w") as f:
                    SeqIO.write(cds_records, f, "fasta")
                fd_a.close()
                fd_c.close()
                aln_input, cds_input = fd_a.name, fd_c.name
            ok, err, failed_id = _backtranslate_pal2nal(aln_input, cds_input, cds_aligned)
            if exclude_for_pal2nal and aln_input != protein_aln:
                try:
                    os.unlink(aln_input)
                    os.unlink(cds_input)
                except OSError:
                    pass
            if ok:
                pal2nal_ok = True
                if exclude_for_pal2nal:
                    from Bio import AlignIO
                    from Bio.Align import MultipleSeqAlignment
                    aln_kept = [r for r in AlignIO.read(protein_aln, "fasta")
                                if r.id not in exclude_for_pal2nal]
                    with open(protein_aln, "w") as f:
                        AlignIO.write(MultipleSeqAlignment(aln_kept), f, "fasta")
                    prot_kept = [r for r in SeqIO.parse(protein_fasta, "fasta")
                                 if r.id not in exclude_for_pal2nal]
                    with open(protein_fasta, "w") as f:
                        SeqIO.write(prot_kept, f, "fasta")
                    cds_kept = [r for r in SeqIO.parse(combined_cds, "fasta")
                                if r.id not in exclude_for_pal2nal]
                    with open(combined_cds, "w") as f:
                        SeqIO.write(cds_kept, f, "fasta")
                    print(
                        f"  {protein}: excluded {len(exclude_for_pal2nal)} sequence(s) "
                        f"with pep/nuc inconsistency."
                    )
                break
            if not failed_id and err:
                failed_id = _parse_pal2nal_failed_id(err)
            if failed_id:
                exclude_for_pal2nal.add(failed_id)
                shown = failed_id[:80] + ("..." if len(failed_id) > 80 else "")
                print(f"  {protein}: excluding pep/nuc inconsistency, retrying without: {shown}")
                continue
            print(f"PAL2NAL failed for {protein}: {(err or 'unknown')[:500]}")
            break
        else:
            print(f"PAL2NAL failed for {protein}: too many inconsistent sequences")
        if not pal2nal_ok:
            continue

        n_in_alignment = len(list(SeqIO.parse(cds_aligned, "fasta")))
        print(f"  {protein}: codon alignment has {n_in_alignment} sequences")

        mafft_prot_dir = os.path.join(alignment_dir, "MAFFT", protein)
        trim_prot_dir = os.path.join(alignment_dir, "Trimmed", protein)
        os.makedirs(mafft_prot_dir, exist_ok=True)
        os.makedirs(trim_prot_dir, exist_ok=True)
        for name, src in [
            ("cds.fasta", combined_cds),
            ("protein.fasta", protein_fasta),
            ("protein_aln.fasta", protein_aln),
            ("cds_aligned.fasta", cds_aligned),
        ]:
            shutil.copy2(src, os.path.join(mafft_prot_dir, name))

        trimmed_fasta = os.path.join(trim_prot_dir, "cds_trimmed.fasta")
        r = subprocess.run(
            ["trimal", "-in", cds_aligned, "-out", trimmed_fasta, "-automated1"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(f"  TrimAl failed for {protein}, copying untrimmed")
            shutil.copy2(cds_aligned, trimmed_fasta)
        else:
            print(f"  Wrote trimmed codon alignment: {trimmed_fasta}")

        # FigTree metadata next to trimmed alignment
        meta_out = os.path.join(trim_prot_dir, f"{protein.lower()}_metadata.txt")
        create_figtree_metadata(
            trimmed_fasta, meta_out,
            segment=LASSA_PROTEIN_REGIONS[protein]['segment'],
        )

    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass
    return extraction_counts


def run_protein_phylogeny(output_dir, proteins, threads=1):
    """Build an IQ-TREE phylogeny for each protein trimmed codon alignment."""
    iqtree_cmd = _resolve_iqtree_cmd()
    if iqtree_cmd is None:
        print("Error: iqtree2/iqtree not found. Required for protein phylogeny.")
        sys.exit(1)

    phylogeny_base = os.path.join(output_dir, "Phylogeny")
    os.makedirs(phylogeny_base, exist_ok=True)
    trimmed_dir = os.path.join(output_dir, "Alignment", "Trimmed")
    mafft_dir = os.path.join(output_dir, "Alignment", "MAFFT")
    nt = "AUTO" if threads <= 0 else str(threads)

    for protein in proteins:
        trimmed_fasta = os.path.join(trimmed_dir, protein, "cds_trimmed.fasta")
        cds_aligned = os.path.join(mafft_dir, protein, "cds_aligned.fasta")
        src = trimmed_fasta if os.path.isfile(trimmed_fasta) else cds_aligned
        if not os.path.isfile(src):
            print(f"  {protein}: no alignment found, skipping phylogeny")
            continue
        protdir = os.path.join(phylogeny_base, protein)
        os.makedirs(protdir, exist_ok=True)
        dest = os.path.join(protdir, "cds_aligned.fasta")
        shutil.copy2(src, dest)
        # Copy metadata if present
        meta_src = os.path.join(trimmed_dir, protein, f"{protein.lower()}_metadata.txt")
        if os.path.isfile(meta_src):
            shutil.copy2(meta_src, os.path.join(protdir, f"{protein.lower()}_metadata.txt"))

        print(f"\nBuilding phylogenetic tree for {protein} using {iqtree_cmd}...")
        cwd = os.getcwd()
        os.chdir(protdir)
        try:
            r = subprocess.run(
                [iqtree_cmd, "-s", "cds_aligned.fasta", "-m", "MFP",
                 "-bb", "10000", "-nt", nt],
                capture_output=True, text=True,
            )
            if r.returncode != 0:
                print(f"IQ-TREE failed for {protein}: {r.stderr}")
            else:
                print(f"IQ-TREE completed for {protein} in Phylogeny/{protein}/")
                for file in os.listdir('.'):
                    if file.endswith(('.bionj', '.ckp.gz', '.mldist', '.log')):
                        try:
                            os.remove(file)
                        except OSError:
                            pass
        finally:
            os.chdir(cwd)


def prompt_alignment_and_proteins(args, interactive):
    """Fill alignment / proteins / phylogeny from CLI or interactive prompts."""
    alignment_from_cli = args.alignment is not None
    phylogeny_from_cli = bool(getattr(args, '_phylogeny_explicit', False))

    if args.alignment is None:
        if interactive:
            print("\nDo you want to perform alignment?")
            print("1. Whole-genome (L and S segment) alignment")
            print("2. Protein (CDS) alignment")
            print("3. No alignment")
            args.alignment = get_user_input(
                "Select alignment option (1-3):", ['1', '2', '3']
            )
        else:
            # Backward compatible: --phylogeny alone implies whole-segment alignment
            args.alignment = '1' if args.phylogeny else '3'

    if args.alignment == '2':
        if not args.proteins:
            if interactive or not alignment_from_cli:
                print("\nSelect proteins for CDS alignment (comma-separated):")
                for num, name in PROTEIN_NUM.items():
                    print(f"  {num} = {name} — {PROTEIN_DESCRIPTIONS[name]}")
                print("  Or use names: L, Z, NP, GPC")
                choice = input("Enter proteins (default: L,NP,GPC,Z): ").strip() or "L,Z,NP,GPC"
                args.proteins = choice
            else:
                args.proteins = "L,Z,NP,GPC"
        proteins = normalize_protein_list(args.proteins)
        invalid = [p for p in proteins if p not in LASSA_PROTEIN_REGIONS]
        if invalid:
            print(f"Error: unknown protein(s): {', '.join(invalid)}")
            print(f"Valid options: {', '.join(LASSA_PROTEIN_REGIONS.keys())} (or numbers 1-4)")
            sys.exit(1)
        args.proteins = ','.join(proteins)

    if args.alignment in ('1', '2') and not phylogeny_from_cli and interactive:
        phy = get_user_input(
            "\nBuild phylogenetic tree(s) from the alignment?\n"
            "1 = Yes\n"
            "2 = No",
            ['1', '2']
        )
        args.phylogeny = (phy == '1')
    elif args.alignment == '1' and args.phylogeny:
        pass  # whole-segment phylogeny as before
    elif args.alignment == '3':
        args.phylogeny = False

    return args


def filter_by_lineage(sequences, target_lineage=None, target_sublineage=None, l_sublineage=None, s_sublineage=None):
    """Filter sequences based on lineage and/or sublineage, with separate L and S sublineage options"""
    if not any([target_lineage, target_sublineage, l_sublineage, s_sublineage]):
        return sequences
        
    filtered_sequences = []
    
    # Get the lineages directory path
    lineages_dir = get_lineages_dir()
    
    # Normalize target lineage to include "Lineage" prefix if not present
    if target_lineage and not target_lineage.startswith("Lineage"):
        target_lineage = f"Lineage{target_lineage}"
    
    # Normalize sublineages to include "Sublineage" prefix if not present
    if target_sublineage and not target_sublineage.startswith("Sublineage"):
        target_sublineage = f"Sublineage{target_sublineage}"
    if l_sublineage and not l_sublineage.startswith("Sublineage"):
        l_sublineage = f"Sublineage{l_sublineage}"
    if s_sublineage and not s_sublineage.startswith("Sublineage"):
        s_sublineage = f"Sublineage{s_sublineage}"
    
    # Load all lineage information first
    l_lineages = {}  # {accession: (lineage, sublineage)}
    s_lineages = {}
    
    # Load L segment lineages
    l_file = os.path.join(lineages_dir, 'l_lineages.txt')
    if os.path.exists(l_file):
        print(f"Reading L segment lineage information from: {l_file}")  # Add debug print
        with open(l_file) as f:
            header = next(f).strip().split('\t')
            acc_idx = next((i for i, col in enumerate(header) if 'acc' in col.lower()), 0)
            lin_idx = next((i for i, col in enumerate(header) if 'lin' in col.lower()), 1)
            sub_idx = next((i for i, col in enumerate(header) if 'sub' in col.lower()), 2)
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    acc = parts[acc_idx].strip()
                    lin = parts[lin_idx].strip()
                    sub = parts[sub_idx].strip()
                    l_lineages[acc] = (lin, sub)
    else:
        print(f"Warning: Could not find L segment lineage file: {l_file}")  # Add debug print
    
    # Load S segment lineages
    s_file = os.path.join(lineages_dir, 's_lineages.txt')
    if os.path.exists(s_file):
        print(f"Reading S segment lineage information from: {s_file}")  # Add debug print
        with open(s_file) as f:
            header = next(f).strip().split('\t')
            acc_idx = next((i for i, col in enumerate(header) if 'acc' in col.lower()), 0)
            lin_idx = next((i for i, col in enumerate(header) if 'lin' in col.lower()), 1)
            sub_idx = next((i for i, col in enumerate(header) if 'sub' in col.lower()), 2)
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    acc = parts[acc_idx].strip()
                    lin = parts[lin_idx].strip()
                    sub = parts[sub_idx].strip()
                    s_lineages[acc] = (lin, sub)
    else:
        print(f"Warning: Could not find S segment lineage file: {s_file}")  # Add debug print
    
    # Now filter sequences using the loaded lineage information
    for seq in sequences:
        record = seq['record']
        accession = record.id.split('_')[0]  # Get accession from the first part
        
        # Handle special cases
        if "Reference" in record.id:
            accession = "NC_004297.1" if 'L_segment' in str(seq['segment']) else "NC_004296.1"
        elif "outgroup" in record.id:
            accession = "KM822127.1" if 'L_segment' in str(seq['segment']) else "KM822128.1"
        
        # Determine if this is an L or S segment
        is_l_segment = 'L_segment' in str(seq['segment']) or seq['segment'] == 'L'
        lineages = l_lineages if is_l_segment else s_lineages
        
        if accession in lineages:
            lineage, sublineage = lineages[accession]
            
            # Check lineage match if specified
            lineage_match = not target_lineage or lineage == target_lineage
            
            # Check sublineage match based on segment type and specified filters
            if is_l_segment:
                sublineage_match = (not l_sublineage and not target_sublineage) or \
                                 (l_sublineage and sublineage == l_sublineage) or \
                                 (not l_sublineage and target_sublineage and sublineage == target_sublineage)
            else:
                sublineage_match = (not s_sublineage and not target_sublineage) or \
                                 (s_sublineage and sublineage == s_sublineage) or \
                                 (not s_sublineage and target_sublineage and sublineage == target_sublineage)
            
            if lineage_match and sublineage_match:
                filtered_sequences.append(seq)
    
    return filtered_sequences

def cli_main():
    try:
        parser = argparse.ArgumentParser(
            description='Download and filter Lassa virus sequences',
            epilog='''example:
  # Download complete genomes from human hosts with known location and date from multiple countries:
  lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --countries "Sierra Leone, Guinea"
  
  # Filter by different sublineages for L and S segments:
  lassaseq -o lassa_output --genome 1 --l_sublineage III --s_sublineage II''')
        
        parser.add_argument('-o', '--outdir', required=True,
                          help='Output directory for sequences')
        
        parser.add_argument('--genome', choices=['1', '2', '3'],
                          help='''Genome completeness filter:
1 = Complete genomes only (>99 percent of reference length)
2 = Partial genomes (specify minimum percent with --completeness)
3 = No completeness filter''')
        
        parser.add_argument('--completeness', type=float,
                          help='''Minimum sequence completeness (1-100 percent)
Required when --genome=2''')
        
        parser.add_argument('--host', choices=['1', '2', '3', '4'],
                          help='''Host filter:
1 = Human sequences only
2 = Rodent sequences only
3 = Both human and rodent sequences
4 = No host filter''')
        
        parser.add_argument('--metadata', choices=['1', '2', '3', '4'],
                          help='''Metadata completeness filter:
1 = Keep only sequences with known location
2 = Keep only sequences with known date
3 = Keep only sequences with both known location and date
4 = No metadata filter''')
        
        parser.add_argument('--countries', type=str,
                          help='''Comma-separated list of countries to filter sequences
Examples: "Sierra Leone, Guinea" or "Nigeria, Mali"
Available: Nigeria, Sierra Leone, Liberia, Guinea, Mali,
          Ghana, Benin, Burkina Faso, Ivory Coast, Togo''')
        
        parser.add_argument('--remove', type=str,
                          help='''(Optional) File containing accession numbers to remove
One accession number per line, lines starting with # are ignored''')
        
        parser.add_argument('--alignment', '-a', type=str, choices=['1', '2', '3'],
                          default=None, metavar='{1,2,3}',
                          help='''Alignment type:
1 = Whole-genome (L and S segment) alignment
2 = Protein (CDS) alignment
3 = No alignment
If omitted with --phylogeny, defaults to 1 (backward compatible)''')

        parser.add_argument('--proteins', '-pr', type=str, default=None, metavar='LIST',
                          help='''For --alignment 2: comma-separated proteins
Names: L,Z,NP,GPC  or numbers: 1=L, 2=Z, 3=NP, 4=GPC
Example: --proteins GPC,NP or -pr 4,3''')

        parser.add_argument('--min-cds-fraction', type=float, default=0.5, metavar='FRAC',
                          help='''For --alignment 2: minimum fraction of RefSeq CDS length
to keep a sequence (default: 0.5)''')

        parser.add_argument('--threads', '-t', type=int, default=1,
                          help='''Threads for minimap2, MAFFT, and IQ-TREE (default: 1).
0 = use all CPUs.''')

        parser.add_argument('--phylogeny', '-p', action='store_true',
                          help='''(Optional) Build phylogenetic tree(s) from the alignment
Whole-genome: one tree per segment. Protein: one tree per protein''')
        
        parser.add_argument('--consensus_L', type=str,
                          help='''(Optional) Path to custom consensus sequence for L segment
The sequence should be in FASTA format''')
        
        parser.add_argument('--consensus_S', type=str,
                          help='''(Optional) Path to custom consensus sequence for S segment
The sequence should be in FASTA format''')
        
        parser.add_argument('--lineage', type=str,
                          help='''(Optional) Filter sequences by lineage
Example: --lineage IV''')
        
        parser.add_argument('--sublineage', type=str,
                          help='''(Optional) Filter sequences by sublineage (applies to both segments)
Example: --sublineage III''')
        
        parser.add_argument('--l_sublineage', type=str,
                          help='''(Optional) Filter L segment sequences by sublineage
Example: --l_sublineage III''')
        
        parser.add_argument('--s_sublineage', type=str,
                          help='''(Optional) Filter S segment sequences by sublineage
Example: --s_sublineage II''')

        # Parse arguments and show help if no arguments provided
        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit(0)
            
        args = parser.parse_args()
        args._phylogeny_explicit = ('--phylogeny' in sys.argv) or ('-p' in sys.argv)
        
        # Initialize countries as None if not provided
        countries = None
        if args.countries:
            countries = [country.strip() for country in args.countries.split(',')]
        
        # Interactive mode if optional arguments are not provided
        interactive = not args.genome
        if not args.genome:
            args.genome = get_user_input(
                "\nSelect genome completeness filter:\n"
                "1 = Complete genomes only (>99 percent)\n"
                "2 = Partial genomes (specify minimum percent)\n"
                "3 = No completeness filter",
                ['1', '2', '3'])
            
            if args.genome == '2' and not args.completeness:
                args.completeness = get_completeness()
        
        if not args.host:
            args.host = get_user_input(
                "\nSelect host filter:\n"
                "1 = Human sequences only\n"
                "2 = Rodent sequences only\n"
                "3 = Both human and rodent sequences\n"
                "4 = No host filter",
                ['1', '2', '3', '4'])
        
        if not args.metadata:
            args.metadata = get_user_input(
                "\nSelect metadata completeness filter:\n"
                "1 = Keep only sequences with known location\n"
                "2 = Keep only sequences with known date\n"
                "3 = Keep only sequences with both known location and date\n"
                "4 = No metadata filter",
                ['1', '2', '3', '4'])

        # Alignment / protein / phylogeny (interactive or CLI defaults)
        args = prompt_alignment_and_proteins(args, interactive)
        
        # Continue with the rest of the code...
        print("\nStarting Lassa virus sequence download")
        print(f"Output directory: {args.outdir}")
        if args.alignment == '2':
            print(f"Alignment: protein CDS ({args.proteins})")
        elif args.alignment == '1':
            print("Alignment: whole-genome (L and S segments)")
        if args.phylogeny:
            print("Phylogeny: enabled")
        
        # Create output directories
        fasta_dir = os.path.join(args.outdir, "FASTA")
        l_segment_dir = os.path.join(fasta_dir, "L_segment")
        s_segment_dir = os.path.join(fasta_dir, "S_segment")
        unknown_segment_dir = os.path.join(fasta_dir, "unknown_segment")
        
        for directory in [fasta_dir, l_segment_dir, s_segment_dir, unknown_segment_dir]:
            os.makedirs(directory, exist_ok=True)

        # Handle consensus sequences if provided
        if args.consensus_L:
            copy_consensus_sequence(args.consensus_L, "L", args.outdir)
        if args.consensus_S:
            copy_consensus_sequence(args.consensus_S, "S", args.outdir)

        # Download reference and outgroup sequences first
        download_and_write_special_sequences(args.outdir)

        # Download sequences
        sequences, segment_counts, location_counts = fetch_sequences()
        initial_total = len(sequences)
        initial_segment_counts = calculate_segment_counts(sequences)

        # After downloading sequences but before filtering, remove specified sequences
        if args.remove and os.path.exists(args.remove):
            with open(args.remove) as f:
                # Read accessions to remove, skip comments and empty lines
                remove_accessions = {line.strip() for line in f 
                                  if line.strip() and not line.startswith('#')}
                
            # Filter out sequences with matching accessions
            original_count = len(sequences)
            sequences = [seq for seq in sequences 
                       if seq['record'].id not in remove_accessions]
            removed_count = original_count - len(sequences)
            
            print(f"\nRemoved {removed_count} sequences based on remove.txt")

        # Apply all filters in sequence
        if args.genome != '3':  # Apply completeness filter if requested
            completeness_filtered = []
            for seq in sequences:
                if args.genome == '1' and is_complete_sequence(seq['record']):
                    completeness_filtered.append(seq)
                elif args.genome == '2' and meets_minimum_completeness(seq['record'], args.completeness):
                    completeness_filtered.append(seq)
        else:
            completeness_filtered = sequences

        # Apply remaining filters
        host_filtered = filter_by_host(completeness_filtered, args.host)
        metadata_filtered = filter_by_metadata(host_filtered, args.metadata)
        
        # Process country filtering
        if countries:
            final_filtered = filter_by_country(metadata_filtered, countries)
        else:
            final_filtered = metadata_filtered

        # After metadata and country filtering, add lineage filtering
        if args.lineage or args.sublineage or args.l_sublineage or args.s_sublineage:
            final_filtered = filter_by_lineage(final_filtered, args.lineage, args.sublineage, 
                                            args.l_sublineage, args.s_sublineage)
            print(f"\nApplied lineage/sublineage filter:")
            if args.lineage:
                print(f"Lineage: {args.lineage}")
            if args.sublineage:
                print(f"Sublineage (both segments): {args.sublineage}")
            if args.l_sublineage:
                print(f"L segment sublineage: {args.l_sublineage}")
            if args.s_sublineage:
                print(f"S segment sublineage: {args.s_sublineage}")
            print(f"Remaining sequences: {len(final_filtered)}")

        # Process final filtered sequences and write FASTA files
        l_sequences, s_sequences, unknown_sequences, unknown_headers = process_sequences(final_filtered, 'both')
        
        # Write main sequence files
        l_output = os.path.join(l_segment_dir, "lassa_l_segments.fasta")
        s_output = os.path.join(s_segment_dir, "lassa_s_segments.fasta")
        unknown_output = os.path.join(unknown_segment_dir, "lassa_unknown_segments.fasta")
        
        SeqIO.write(l_sequences, l_output, "fasta")
        SeqIO.write(s_sequences, s_output, "fasta")
        SeqIO.write(unknown_sequences, unknown_output, "fasta")
        
        # Update the output structure in summary
        written_counts = {
            'L': len(l_sequences),
            'S': len(s_sequences),
            'unknown': len(unknown_sequences)
        }

        # Write summary file
        write_summary(args.outdir, initial_total, len(host_filtered), 
                     initial_segment_counts, location_counts, 'both', 
                     written_counts, host_filtered, args.genome, 
                     args.completeness, sequences, args.host, completeness_filtered)
        
        # Append metadata and country filtering results
        with open(os.path.join(args.outdir, 'summary_Lassa.txt'), 'a') as f:
            write_metadata_filtering_summary(f, host_filtered, metadata_filtered, args.metadata, written_counts)
            write_country_filtering_summary(f, metadata_filtered, final_filtered, countries, written_counts)

        print(f"\nFinal counts after all filtering:")
        print(f"Wrote {written_counts['L']} L segments, {written_counts['S']} S segments, "
              f"and {written_counts['unknown']} unknown segments")
        print(f"Updated summary report written to: {os.path.join(args.outdir, 'summary_Lassa.txt')}")

        # Alignment and optional phylogeny
        nthreads = args.threads if getattr(args, 'threads', 1) > 0 else (os.cpu_count() or 1)

        if args.alignment == '1':
            print("\nCreating files for whole-genome (segment) alignment...")
            phylogeny_dir = create_phylogeny_directories(args.outdir)
            l_count = concatenate_fasta_files(args.outdir, phylogeny_dir, "L")
            s_count = concatenate_fasta_files(args.outdir, phylogeny_dir, "S")
            print(f"\nCreated concatenated FASTA files:")
            print(f"L segment: {l_count} unique sequences")
            print(f"S segment: {s_count} unique sequences")
            perform_phylogenetic_analysis(phylogeny_dir, "L", do_tree=args.phylogeny)
            perform_phylogenetic_analysis(phylogeny_dir, "S", do_tree=args.phylogeny)
            print("\nWhole-genome analysis completed. Output files:")
            print(f"Alignments: {os.path.join(phylogeny_dir, 'MSA')}")
            print(f"Trimmed alignments: {os.path.join(phylogeny_dir, 'TrimAl')}")
            if args.phylogeny:
                print(f"Tree files: {os.path.join(phylogeny_dir, 'Tree')}")

        elif args.alignment == '2':
            proteins = normalize_protein_list(args.proteins)
            print(f"\nRunning protein (CDS) alignment for: {', '.join(proteins)}")
            run_protein_cds_pipeline(
                args.outdir,
                proteins,
                min_cds_fraction=getattr(args, 'min_cds_fraction', 0.5),
                threads=nthreads,
            )
            print(f"\nProtein alignment completed. Output: {os.path.join(args.outdir, 'Alignment')}")
            if args.phylogeny:
                print("\nBuilding protein phylogenies...")
                run_protein_phylogeny(args.outdir, proteins, threads=nthreads)
                print(f"Protein phylogenies completed. Output: {os.path.join(args.outdir, 'Phylogeny')}")

    except Exception as e:
        print(f"\nAn unexpected error occurred: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    cli_main() 
