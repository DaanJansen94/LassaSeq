# LassaSeq

A simple tool to download and organize Lassa virus sequences from GenBank.

## Installation

```bash
   pip install lassaseq
```

## Usage

```bash
lassaseq -o output_directory
```

### Options

```bash
- `-o, --outdir`: Output directory for sequences (required)
- `--genome`: Genome completeness filter
  - 1: Complete genomes only (>99%)
  - 2: Partial genomes (specify minimum completeness)
  - 3: No completeness filter
- `--completeness`: Minimum sequence completeness (1-100), required when --genome=2
- `--host`: Host filter
  - 1: Human sequences
  - 2: Non-human sequences
  - 3: No host filter
```

### Output

- Organized FASTA files in separate directories for L and S segments
- Summary file with:
  - Host distribution
  - Segment counts
  - Geographical distribution
  - Filtering results

## Example

```bash
lassaseq -o lassa_sequences --genome 1 --host 1 
lassaseq -o lassa_sequences --genome 2 --completeness 80 --host 3 
```

## Requirements
- Python 3.6+
- Biopython