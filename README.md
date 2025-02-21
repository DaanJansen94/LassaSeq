# LassaSeq
LassaSeq is a command-line tool that simplifies the process of analyzing Lassa virus sequences. It automates the complete workflow from downloading sequences to creating phylogenetic trees, with special handling for Lassa's bi-segmented genome.

## Genome Segments
Lassa virus has a bi-segmented RNA genome consisting of:
- **L segment**: (~7.2kb) Encodes the RNA-dependent RNA polymerase and Z protein
- **S segment**: (~3.4kb) Encodes the nucleoprotein (NP) and glycoprotein precursor (GPC)


## Installation

1. First, install conda if you haven't already:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Create and activate a new conda environment:
   ```bash
   conda create -n lassa_env -c bioconda python=3.9 mafft trimal iqtree
   conda activate lassa_env
   ```

3. Install LassaSeq:
   ```bash
   git clone https://github.com/DaanJansen94/LassaSeq.git   
   cd LassaSeq
   pip install .
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
