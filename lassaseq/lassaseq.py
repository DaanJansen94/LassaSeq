from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm
import argparse
import os
from datetime import datetime
import calendar
import time
import sys

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

def clean_country_name(country):
    """Clean country name to remove spaces and special characters"""
    # Dictionary for specific country name fixes
    country_fixes = {
        'Sierra Leone': 'SierraLeone',
        'Burkina Faso': 'BurkinaFaso',
        'Costa Rica': 'CostaRica',
        'South Africa': 'SouthAfrica',
        "Cote d'Ivoire": 'IvoryCoast',
        "Côte d'Ivoire": 'IvoryCoast',
        'Ivory Coast': 'IvoryCoast'
        # Add more multi-word countries as needed
    }
    
    # Check if country is in our fixes dictionary
    if country in country_fixes:
        return country_fixes[country]
    
    # Otherwise, remove spaces and special characters
    return country.replace(' ', '').replace('-', '').replace(',', '')

def get_metadata(record):
    """Extract metadata from record"""
    # Get accession
    accession = record.id
    
    # Initialize location and date
    location = "UnknownLoc"
    collection_date = "UnknownDate"
    
    # Look for location and date in source features
    for feature in record.features:
        if feature.type == "source":
            # Get location from geo_loc_name
            if 'geo_loc_name' in feature.qualifiers:
                geo_loc = feature.qualifiers['geo_loc_name'][0]
                # Only use location if it's not 'missing'
                if 'missing' not in geo_loc.lower():
                    # Extract country (part before the colon)
                    location = geo_loc.split(':')[0].strip()
                    # Clean up location string
                    location = clean_country_name(location)
            
            # Get collection date
            if 'collection_date' in feature.qualifiers:
                date_str = feature.qualifiers['collection_date'][0]
                if date_str.lower() != 'missing':
                    collection_date = convert_date_to_decimal_year(date_str)
    
    return f"{accession}_{location}_{collection_date}"
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

def write_summary(outdir, initial_total, filtered_total, initial_segment_counts, location_counts, requested_segment, written_counts, filtered_sequences, genome_choice, completeness=None, sequences=None, host_choice=None):
    """Write summary report to a file"""
    summary_file = os.path.join(outdir, "summary_Lassa.txt")
    
    with open(summary_file, 'w') as f:
        f.write("Lassa Virus Sequence Download Summary\n")
        f.write("====================================\n\n")
        f.write(f"Total Lassa virus sequences found: {initial_total}\n")
        
        # Initial distribution
        f.write("\nInitial Segment Distribution:\n")
        f.write("---------------------------\n")
        f.write(f"L segments: {initial_segment_counts['L']}\n")
        f.write(f"S segments: {initial_segment_counts['S']}\n")
        f.write(f"Unknown segments: {initial_segment_counts['unknown']}\n\n")
        
        f.write("Initial Geographical Distribution:\n")
        f.write("--------------------------------\n")
        f.write("Country        L segments    S segments    Unknown         Total\n")
        f.write("-----------   -----------   -----------   ---------   ----------\n")
        
        # Calculate initial location counts
        initial_locations = {'L': {}, 'S': {}, 'unknown': {}}
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
                            location = clean_country_name(location)
            
            if location not in initial_locations[segment_type]:
                initial_locations[segment_type][location] = 0
            initial_locations[segment_type][location] += 1
        
        write_geographical_distribution(f, initial_locations)
        
        # After completeness filtering (if applied)
        if genome_choice != '3':
            f.write("\nAfter Completeness Filtering:\n")
            f.write("--------------------------\n")
            if genome_choice == '1':
                f.write(f"Filtering for completeness (>99%)\n")
            else:
                f.write(f"Filtering for completeness (>{completeness}%)\n")
            
            completeness_locations = {'L': {}, 'S': {}, 'unknown': {}}
            completeness_counts = {'L': 0, 'S': 0, 'unknown': 0}
            
            for seq in filtered_sequences:
                record = seq['record']
                segment = get_segment_type(record)
                segment_type = 'unknown' if segment is None else segment
                completeness_counts[segment_type] += 1
                
                location = "UnknownLoc"
                for feature in record.features:
                    if feature.type == "source":
                        if 'geo_loc_name' in feature.qualifiers:
                            geo_loc = feature.qualifiers['geo_loc_name'][0]
                            if 'missing' not in geo_loc.lower():
                                location = geo_loc.split(':')[0].strip()
                                location = clean_country_name(location)
                
                if location not in completeness_locations[segment_type]:
                    completeness_locations[segment_type][location] = 0
                completeness_locations[segment_type][location] += 1
            
            f.write(f"\nSegment counts after completeness filtering:\n")
            f.write(f"L segments: {completeness_counts['L']}\n")
            f.write(f"S segments: {completeness_counts['S']}\n")
            f.write(f"Unknown segments: {completeness_counts['unknown']}\n\n")
            
            f.write("Geographical Distribution After Completeness Filtering:\n")
            f.write("------------------------------------------------\n")
            f.write("Country        L segments    S segments    Unknown         Total\n")
            f.write("-----------   -----------   -----------   ---------   ----------\n")
            write_geographical_distribution(f, completeness_locations)
        
        # After host filtering (if applied)
        if host_choice != '3':
            f.write("\nAfter Host Filtering:\n")
            f.write("-------------------\n")
            if host_choice == '1':
                f.write("Filtering for human host only\n")
            else:
                f.write("Filtering for non-human host only\n")
            
            final_locations = {'L': {}, 'S': {}, 'unknown': {}}
            final_counts = {'L': 0, 'S': 0, 'unknown': 0}
            
            for seq in filtered_sequences:
                record = seq['record']
                segment = get_segment_type(record)
                segment_type = 'unknown' if segment is None else segment
                final_counts[segment_type] += 1
                
                location = "UnknownLoc"
                for feature in record.features:
                    if feature.type == "source":
                        if 'geo_loc_name' in feature.qualifiers:
                            geo_loc = feature.qualifiers['geo_loc_name'][0]
                            if 'missing' not in geo_loc.lower():
                                location = geo_loc.split(':')[0].strip()
                                location = clean_country_name(location)
                
                if location not in final_locations[segment_type]:
                    final_locations[segment_type][location] = 0
                final_locations[segment_type][location] += 1
            
            f.write(f"\nFinal segment counts:\n")
            f.write(f"L segments: {final_counts['L']}\n")
            f.write(f"S segments: {final_counts['S']}\n")
            f.write(f"Unknown segments: {final_counts['unknown']}\n\n")
            
            f.write("Final Geographical Distribution:\n")
            f.write("-----------------------------\n")
            f.write("Country        L segments    S segments    Unknown         Total\n")
            f.write("-----------   -----------   -----------   ---------   ----------\n")
            write_geographical_distribution(f, final_locations)
        
        # Write output files section
        f.write("\nOutput files:\n")
        f.write("  FASTA/L_segment/lassa_l_segments.fasta\n")
        f.write("  FASTA/S_segment/lassa_s_segments.fasta\n")
        f.write("  FASTA/unknown_segment/lassa_unknown_segments.fasta\n")

def write_geographical_distribution(f, locations):
    """Helper function to write geographical distribution tables"""
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

def cli_main():
    try:
        parser = argparse.ArgumentParser(description='Download Lassa virus sequences')
        parser.add_argument('-o', '--outdir', required=True, help='Output directory for sequences')
        
        # Modified genome completeness options
        parser.add_argument('--genome', choices=['1', '2', '3'],
                           help='Genome completeness: 1=Complete only (>99%), 2=Partial (requires --completeness), 3=No filter')
        parser.add_argument('--completeness', type=float,
                           help='Minimum sequence completeness (1-100), required when --genome=2')
        parser.add_argument('--host', choices=['1', '2', '3'],
                           help='Host filter: 1=Human only, 2=Non-human only, 3=No host filter')
        
        args = parser.parse_args()
        
        print("\nStarting Lassa virus sequence download")
        print(f"Output directory: {args.outdir}")
        
        # Create segment-specific directories
        fasta_dir = os.path.join(args.outdir, "FASTA")
        l_segment_dir = os.path.join(fasta_dir, "L_segment")
        s_segment_dir = os.path.join(fasta_dir, "S_segment")
        unknown_segment_dir = os.path.join(fasta_dir, "unknown_segment")
        
        # Create all directories
        for directory in [fasta_dir, l_segment_dir, s_segment_dir, unknown_segment_dir]:
            os.makedirs(directory, exist_ok=True)
        
        # Handle genome completeness options
        if args.genome:
            genome_choice = args.genome
            completeness = args.completeness if args.genome == '2' else None
        else:
            print("\nGenome completeness options:")
            print("1. Complete genomes only (>99%)")
            print("2. Partial genomes (specify minimum completeness)")
            print("3. No completeness filter")
            genome_choice = get_user_input("", ['1', '2', '3'])
            
            completeness = None
            if genome_choice == '2':
                completeness = get_completeness()
        
        # Handle host filtering options
        if args.host:
            host_choice = args.host
        else:
            print("\nHost filtering options:")
            print("1. Human sequences only")
            print("2. Non-human sequences only")
            print("3. No host filter")
            host_choice = get_user_input("", ['1', '2', '3'])
        
        sequences, segment_counts, location_counts = fetch_sequences()
        initial_total = len(sequences)
        
        # Get initial counts before filtering
        initial_segment_counts = {'L': 0, 'S': 0, 'unknown': 0}
        for seq in sequences:
            segment_type = get_segment_type(seq['record'])
            if segment_type:
                initial_segment_counts[segment_type] += 1
            else:
                initial_segment_counts['unknown'] += 1
        
        # Filter sequences based on completeness and host
        filtered_sequences = []
        for seq in sequences:
            include_sequence = True
            record = seq['record']
            
            # Completeness filter
            if genome_choice == '1':
                if not is_complete_sequence(record):
                    include_sequence = False
            elif genome_choice == '2':
                if not meets_minimum_completeness(record, completeness):
                    include_sequence = False
            
            # Host filter
            if include_sequence and host_choice != '3':
                host_found = False
                is_human = False
                for feature in record.features:
                    if feature.type == "source":
                        if 'host' in feature.qualifiers:
                            host_found = True
                            host = feature.qualifiers['host'][0].lower()
                            is_human = 'homo sapiens' in host or 'human' in host
                            break
                
                if host_choice == '1':  # Human only
                    if not (host_found and is_human):
                        include_sequence = False
                elif host_choice == '2':  # Non-human only
                    if not (host_found and not is_human):
                        include_sequence = False
            
            if include_sequence:
                filtered_sequences.append(seq)
        
        # Process both segments
        l_sequences, s_sequences, unknown_sequences, unknown_headers = process_sequences(filtered_sequences, 'both')
        
        # Write sequences to their respective directories
        l_output = os.path.join(l_segment_dir, "lassa_l_segments.fasta")
        SeqIO.write(l_sequences, l_output, "fasta")
        
        s_output = os.path.join(s_segment_dir, "lassa_s_segments.fasta")
        SeqIO.write(s_sequences, s_output, "fasta")
        
        unknown_output = os.path.join(unknown_segment_dir, "lassa_unknown_segments.fasta")
        SeqIO.write(unknown_sequences, unknown_output, "fasta")
        
        written_counts = {
            'L': len(l_sequences),
            'S': len(s_sequences),
            'unknown': len(unknown_sequences)
        }
        
        print(f"\nFound {initial_segment_counts['L']} L segments, {initial_segment_counts['S']} S segments, "
              f"and {initial_segment_counts['unknown']} unknown segments before filtering")
        print(f"Wrote {written_counts['L']} L segments, {written_counts['S']} S segments, "
              f"and {written_counts['unknown']} unknown segments after filtering")
        
        # Update the summary file to reflect new directory structure
        write_summary(args.outdir, initial_total, len(filtered_sequences), initial_segment_counts, 
                     location_counts, 'both', written_counts, filtered_sequences, 
                     genome_choice, completeness, sequences, host_choice)
        print(f"Summary report written to: {os.path.join(args.outdir, 'summary_Lassa.txt')}")
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    cli_main() 