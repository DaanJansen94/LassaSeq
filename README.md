# LassaSeq

A command-line tool for downloading and analyzing Lassa virus sequences.

## Installation
pip install git+https://github.com/yourusername/lassaseq.git

## Requirements

- Python 3.6 or higher
- Biopython
- MAFFT (for alignment)
- TrimAl (for alignment trimming)
- IQ-TREE 2 (for phylogenetic analysis)

## Usage

lassaseq --output-dir results --segment L --genome 1 --host 1 --metadata 3 --beast 2 --phylogeny

### Required Arguments

- `--output-dir`: Output directory for results
- `--segment`: Viral segment to analyze (L or S)

### Optional Arguments

- `--genome`: Genome completeness
  - 1 = Complete genomes only
  - 2 = Partial genomes (requires --completeness)
  - 3 = All genomes
- `--completeness`: Required when --genome=2 (value between 1-100)
- `--host`: Host filter
  - 1 = Human only
  - 2 = Non-human only
  - 3 = All hosts
- `--metadata`: Metadata filter
  - 1 = Location only
  - 2 = Date only
  - 3 = Both location and date
  - 4 = None
- `--beast`: Required when --metadata is 2 or 3
  - 1 = No
  - 2 = Yes
- `--consensus-file`: Path to consensus FASTA file to include
- `--phylogeny`: Create phylogenetic tree using IQ-TREE 2
- `--remove`: Path to text file containing headers/accession IDs to remove

## License

MIT License
