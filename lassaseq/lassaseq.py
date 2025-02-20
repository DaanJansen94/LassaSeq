from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm
import argparse
import os
from datetime import datetime
import calendar
import time

def get_segment_type(record):
    """Helper function to determine segment type from various GenBank fields"""
    # First try official segment annotation
    for feature in record.features:
        if feature.type == "source":
            if 'segment' in feature.qualifiers:
                segment_value = feature.qualifiers['segment'][0].upper()
                if segment_value in ['L', 'S']:
                    return segment_value
    
    # Check description and features for protein information
    description = record.description.lower()
    
    # Expanded protein lists
    s_proteins = [
        'nucleoprotein', 'glycoprotein', 'gpc', 'np', 'nuc', 
        'nucleocapsid', 'g protein', 'g1', 'g2', 'gp1', 'gp2'
    ]
    l_proteins = [
        'polymerase', 'z protein', 'l protein', 'rna-dependent', 
        'rna dependent', 'rdrp', 'zinc-binding', 'zinc binding', 
        'matrix protein', 'ring protein'
    ]
    
    # Check description
    for protein in s_proteins:
        if protein in description:
            return 'S'
    for protein in l_proteins:
        if protein in description:
            return 'L'
    
    # Check all features for protein information
    for feature in record.features:
        if feature.type in ['CDS', 'mat_peptide', 'gene']:
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0].lower()
                for protein in s_proteins:
                    if protein in product:
                        return 'S'
                for protein in l_proteins:
                    if protein in product:
                        return 'L'
            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers['gene'][0].lower()
                for protein in s_proteins:
                    if protein in gene:
                        return 'S'
                for protein in l_proteins:
                    if protein in gene:
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
                    segment_type = get_segment_type(record)
                    # Store original header
                    original_header = f"{record.description}"
                    # Add sequence regardless of segment type
                    sequences.append({
                        'id': record.id,
                        'record': record,
                        'segment': segment_type,  # This will be None for unknown segments
                        'original_header': original_header  # Store original header
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

def write_summary(outdir, initial_total, filtered_total, segment_counts, location_counts, requested_segment, written_counts, sequences, genome_choice, completeness=None):
    """Write summary report to a file"""
    summary_file = os.path.join(outdir, "summary_Lassa.txt")
    
    with open(summary_file, 'w') as f:
        f.write("Lassa Virus Sequence Download Summary\n")
        f.write("====================================\n\n")
        f.write(f"Total Lassa virus sequences found: {initial_total}\n")
        
        # Make filtering text dynamic based on choice
        if genome_choice == '1':
            f.write(f"Sequences after filtering for completeness (>99%): {filtered_total}\n\n")
        elif genome_choice == '2':
            f.write(f"Sequences after filtering for completeness (>{completeness}%): {filtered_total}\n\n")
        else:
            f.write(f"No completeness filtering applied: {filtered_total}\n\n")
        
        # Calculate total segments found after filtering
        segment_counts = {'L': 0, 'S': 0, 'unknown': 0}
        for seq in sequences:
            segment = get_segment_type(seq['record'])
            if segment == 'L':
                segment_counts['L'] += 1
            elif segment == 'S':
                segment_counts['S'] += 1
            else:
                segment_counts['unknown'] += 1
        
        f.write(f"Segments found after filtering:\n")
        f.write(f"  L segments: {segment_counts['L']}\n")
        f.write(f"  S segments: {segment_counts['S']}\n")
        f.write(f"  Unknown/unspecified: {segment_counts['unknown']}\n\n")
        
        f.write("Geographical Distribution of Segments:\n")
        f.write("-------------------------------------\n")
        f.write("Country        L segments    S segments    Unknown         Total\n")
        f.write("-----------   -----------   -----------   ---------   ----------\n")
        
        # Recalculate location counts after reclassification
        location_counts = {'L': {}, 'S': {}, 'unknown': {}}
        
        for seq in sequences:
            record = seq['record']
            # Get location
            location = "UnknownLoc"
            for feature in record.features:
                if feature.type == "source":
                    if 'geo_loc_name' in feature.qualifiers:
                        geo_loc = feature.qualifiers['geo_loc_name'][0]
                        if 'missing' not in geo_loc.lower():
                            location = geo_loc.split(':')[0].strip()
                            location = clean_country_name(location)
            
            # Get segment type
            segment = get_segment_type(record)
            
            # Update counts
            if segment == 'L':
                if location not in location_counts['L']:
                    location_counts['L'][location] = 0
                location_counts['L'][location] += 1
            elif segment == 'S':
                if location not in location_counts['S']:
                    location_counts['S'][location] = 0
                location_counts['S'][location] += 1
            else:
                if location not in location_counts['unknown']:
                    location_counts['unknown'][location] = 0
                location_counts['unknown'][location] += 1
        
        # Get all unique countries
        all_countries = set()
        for segment_type in ['L', 'S', 'unknown']:
            all_countries.update(location_counts[segment_type].keys())
        
        # Calculate and write counts for each country
        total_l = 0
        total_s = 0
        total_unknown = 0
        for country in sorted(all_countries):
            l_count = location_counts['L'].get(country, 0)
            s_count = location_counts['S'].get(country, 0)
            unknown_count = location_counts['unknown'].get(country, 0)
            country_total = l_count + s_count + unknown_count
            
            if country_total > 0:  # Only write countries with sequences
                f.write(f"{country:<13} {l_count:>11}   {s_count:>11}   {unknown_count:>9}   {country_total:>10}\n")
            
            total_l += l_count
            total_s += s_count
            total_unknown += unknown_count
        
        # Write totals
        f.write("-" * 60 + "\n")
        total_all = total_l + total_s + total_unknown
        f.write(f"{'Total':<13} {total_l:>11}   {total_s:>11}   {total_unknown:>9}   {total_all:>10}\n\n")
        
        f.write("Sequences downloaded:\n")
        if requested_segment.upper() == 'BOTH':
            f.write(f"  L segments written: {written_counts['L']}\n")
            f.write(f"  S segments written: {written_counts['S']}\n")
            f.write(f"  Unknown segments written: {written_counts['unknown']}\n")
            f.write(f"Output files:\n")
            f.write(f"  L segments: lassa_l_segments.fasta\n")
            f.write(f"  S segments: lassa_s_segments.fasta\n")
            f.write(f"  Unknown segments: lassa_unknown_segments.fasta\n")
        else:
            # Use the correct key (L or S) based on requested segment
            f.write(f"  {requested_segment.upper()} segments written: {written_counts[requested_segment.upper()]}\n")
            f.write(f"Output files:\n")
            f.write(f"  {requested_segment.upper()} segments: lassa_{requested_segment.lower()}_segments.fasta\n")

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
    parser = argparse.ArgumentParser(description='Download Lassa virus sequences')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory for sequences')
    parser.add_argument('-s', '--segment', required=True, choices=['L', 'S', 'both'], 
                       help='Segment type to download (L, S, or both)')
    
    # Modified genome completeness options
    parser.add_argument('--genome', choices=['1', '2'],
                       help='Genome completeness: 1=Complete only (>99%), 2=Partial (requires --completeness)')
    parser.add_argument('--completeness', type=float,
                       help='Minimum sequence completeness (1-100), required when --genome=2')
    
    args = parser.parse_args()
    
    print(f"\nStarting Lassa virus sequence download for {args.segment} segment(s)")
    print(f"Output directory: {args.outdir}")
    
    # Create FASTA directory
    fasta_dir = os.path.join(args.outdir, "FASTA")
    os.makedirs(fasta_dir, exist_ok=True)
    
    # Handle genome completeness options
    if args.genome:
        # Non-interactive mode
        if args.genome == '2' and args.completeness is None:
            parser.error("--completeness is required when --genome=2")
        if args.genome == '2' and not (1 <= args.completeness <= 100):
            parser.error("--completeness must be between 1 and 100")
        genome_choice = args.genome
        completeness = args.completeness if args.genome == '2' else None
    else:
        # Interactive mode
        print("\nGenome completeness options:")
        print("1. Complete genomes only (>99%)")
        print("2. Partial genomes (specify minimum completeness)")
        genome_choice = get_user_input("", ['1', '2'])
        
        completeness = None
        if genome_choice == '2':
            completeness = get_completeness()
    
    sequences, segment_counts, location_counts = fetch_sequences()
    
    # Store initial count
    initial_total = len(sequences)
    
    # Filter sequences based on completeness criteria
    filtered_sequences = []
    for seq in sequences:
        include_sequence = True
        if genome_choice == '1':
            if not is_complete_sequence(seq['record']):
                include_sequence = False
        elif genome_choice == '2':
            if not meets_minimum_completeness(seq['record'], completeness):
                include_sequence = False
        
        if include_sequence:
            filtered_sequences.append(seq)
    
    # Update sequences list with filtered sequences
    sequences = filtered_sequences
    
    if args.segment.upper() == 'BOTH':
        l_sequences, s_sequences, unknown_sequences, unknown_headers = process_sequences(sequences, 'both')
        
        # Write L segments
        l_output = os.path.join(fasta_dir, "lassa_l_segments.fasta")
        SeqIO.write(l_sequences, l_output, "fasta")
        
        # Write S segments
        s_output = os.path.join(fasta_dir, "lassa_s_segments.fasta")
        SeqIO.write(s_sequences, s_output, "fasta")
        
        # Write unknown segments
        unknown_output = os.path.join(fasta_dir, "lassa_unknown_segments.fasta")
        SeqIO.write(unknown_sequences, unknown_output, "fasta")
        
        written_counts = {
            'L': len(l_sequences), 
            'S': len(s_sequences),
            'unknown': len(unknown_sequences)
        }
        print(f"\nFound {len(l_sequences)} L segments, {len(s_sequences)} S segments, and {len(unknown_sequences)} unknown segments before filtering")
        print(f"Wrote {written_counts['L']} L segments, {written_counts['S']} S segments, and {written_counts['unknown']} unknown segments after filtering")
    else:
        filtered_sequences = process_sequences(sequences, args.segment)
        output_file = os.path.join(fasta_dir, f"lassa_{args.segment.lower()}_segments.fasta")
        SeqIO.write(filtered_sequences, output_file, "fasta")
        
        written_counts = {
            'L': len(filtered_sequences) if args.segment.upper() == 'L' else 0,
            'S': len(filtered_sequences) if args.segment.upper() == 'S' else 0,
            'unknown': 0
        }
        segment_type = args.segment.upper()
        total_before = sum(1 for seq in sequences if get_segment_type(seq['record']) == segment_type)
        print(f"\nFound {total_before} {segment_type} segments before filtering")
        print(f"Wrote {len(filtered_sequences)} {segment_type} segments after filtering")
    
    # Write summary report
    write_summary(args.outdir, initial_total, len(sequences), segment_counts, location_counts, 
                 args.segment, written_counts, sequences, genome_choice, completeness)
    print(f"Summary report written to: {os.path.join(args.outdir, 'summary_Lassa.txt')}")

if __name__ == "__main__":
    cli_main() 