# LassaSeq

A simple tool to download and organize Lassa virus sequences from GenBank.

## Installation

```bash
   pip install lassaseq
```

## Usage

Download Lassa virus sequences by segment type:

```bash
lassaseq -o /path/to/output -s L # Download L segments
lassaseq -o /path/to/output -s S # Download S segments
lassaseq -o /path/to/output -s both # Download both L and S segments
```

## Output

- FASTA files with standardized headers: `accession_country_date`
- Summary report with:
  - Total sequences found
  - Segment counts (L/S)
  - Geographical distribution table
  - Download statistics

## Requirements
- Python 3.6+
- Biopython
