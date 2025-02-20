from Bio import Entrez
from Bio import SeqIO
from tqdm import tqdm
import argparse
import os
from datetime import datetime
import calendar

def get_segment_type(record):
    """Helper function to determine segment type from various GenBank fields"""
    for feature in record.features:
        if feature.type == "source":
            if 'segment' in feature.qualifiers:
                segment_value = feature.qualifiers['segment'][0].upper()
                if segment_value in ['L', 'S']:
                    return segment_value
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
    
    print("Searching for Lassa virus sequences...")
    search_term = "Mammarenavirus lassaense[Organism]"
    
    handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    
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
            
        except Exception as e:
            print(f"Error downloading batch starting at {start}: {str(e)}")
            continue
    
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

def write_summary(outdir, total_count, segment_counts, location_counts, requested_segment, written_counts, sequences):
    """Write summary report to a file"""
    summary_file = os.path.join(outdir, "lassa_download_summary.txt")
    
    # Recalculate location counts after reclassification
    final_location_counts = {'L': {}, 'S': {}}
    unknown_locations = {}
    
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
        
        # First try official annotation, then try header-based classification
        final_segment = None
        if seq['segment'] in ['L', 'S']:
            final_segment = seq['segment']
        else:
            # For unknown segments, try to determine from protein names using original header
            final_segment = get_segment_from_protein(seq['original_header'])
        
        # Update counts based on final classification
        if final_segment == 'L':
            if location not in final_location_counts['L']:
                final_location_counts['L'][location] = 0
            final_location_counts['L'][location] += 1
        elif final_segment == 'S':
            if location not in final_location_counts['S']:
                final_location_counts['S'][location] = 0
            final_location_counts['S'][location] += 1
        else:  # Only truly unknown sequences
            if location not in unknown_locations:
                unknown_locations[location] = 0
            unknown_locations[location] += 1
    
    with open(summary_file, 'w') as f:
        f.write("Lassa Virus Sequence Download Summary\n")
        f.write("====================================\n\n")
        f.write(f"Total Lassa virus sequences found: {total_count}\n")
        f.write(f"Segments found in database:\n")
        f.write(f"  L segments: {written_counts['L']}\n")
        f.write(f"  S segments: {written_counts['S']}\n")
        f.write(f"  Unknown/unspecified: {written_counts['unknown']}\n\n")
        
        f.write("Geographical Distribution of Segments:\n")
        f.write("-------------------------------------\n")
        f.write("Country        L segments    S segments    Unknown         Total\n")
        f.write("-----------   -----------   -----------   ---------   ----------\n")
        
        # Get all unique countries
        all_countries = set(final_location_counts['L'].keys()) | \
                       set(final_location_counts['S'].keys()) | \
                       set(unknown_locations.keys())
        
        # Calculate and write counts for each country
        total_l = 0
        total_s = 0
        total_unknown = 0
        for country in sorted(all_countries):
            l_count = final_location_counts['L'].get(country, 0)
            s_count = final_location_counts['S'].get(country, 0)
            unknown_count = unknown_locations.get(country, 0)
            country_total = l_count + s_count + unknown_count
            if country_total > 0:  # Only write countries with sequences
                f.write(f"{country:<13} {l_count:>11}   {s_count:>11}   {unknown_count:>9}   {country_total:>10}\n")
            total_l += l_count
            total_s += s_count
            total_unknown += unknown_count
        
        # Write totals with matching format
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

def cli_main():
    parser = argparse.ArgumentParser(description='Download Lassa virus sequences')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory for sequences')
    parser.add_argument('-s', '--segment', required=True, choices=['L', 'S', 'both'], 
                       help='Segment type to download (L, S, or both)')
    args = parser.parse_args()
    
    print(f"\nStarting Lassa virus sequence download for {args.segment} segment(s)")
    print(f"Output directory: {args.outdir}")
    
    os.makedirs(args.outdir, exist_ok=True)
    
    sequences, segment_counts, location_counts = fetch_sequences()
    
    # Always process both segments for consistent counting
    l_sequences, s_sequences, unknown_sequences, unknown_headers = process_sequences(sequences, 'both')
    
    written_counts = {
        'L': len(l_sequences), 
        'S': len(s_sequences),
        'unknown': len(unknown_sequences)
    }
    
    # Write only the requested segment files
    if args.segment.upper() == 'BOTH':
        # Write L segments
        l_output = os.path.join(args.outdir, "lassa_l_segments.fasta")
        SeqIO.write(l_sequences, l_output, "fasta")
        
        # Write S segments
        s_output = os.path.join(args.outdir, "lassa_s_segments.fasta")
        SeqIO.write(s_sequences, s_output, "fasta")
        
        # Write unknown segments (only for 'both' option)
        unknown_output = os.path.join(args.outdir, "lassa_unknown_segments.fasta")
        SeqIO.write(unknown_sequences, unknown_output, "fasta")
        
        print(f"\nWrote {len(l_sequences)} L segments, {len(s_sequences)} S segments, and {len(unknown_sequences)} unknown segments")
    else:
        # Write only the requested segment
        if args.segment.upper() == 'L':
            output_file = os.path.join(args.outdir, "lassa_l_segments.fasta")
            SeqIO.write(l_sequences, output_file, "fasta")
            print(f"\nWrote {len(l_sequences)} L segments")
        else:  # S segment
            output_file = os.path.join(args.outdir, "lassa_s_segments.fasta")
            SeqIO.write(s_sequences, output_file, "fasta")
            print(f"\nWrote {len(s_sequences)} S segments")
    
    # Write summary report
    write_summary(args.outdir, len(sequences), segment_counts, location_counts, args.segment, written_counts, sequences)
    print(f"Summary report written to: {os.path.join(args.outdir, 'lassa_download_summary.txt')}")

if __name__ == "__main__":
    cli_main() 