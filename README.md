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

## Re-installation
When updates are pushed to GitHub, or when you want to use your own modifications to the code, you'll need to reinstall the package:

```bash
conda activate lassa_env  # Make sure you're in the right environment
cd LassaSeq
git pull  # Get the latest updates from GitHub
pip uninstall LassaSeq
pip install .
```

Note: Any time you modify the code or pull updates from GitHub, you need to reinstall the package using these commands for the changes to take effect.

## Usage

To see all available options:
```bash
lassaseq --help
```

This will display:
```
usage: lassaseq [-h] -o  [--genome {1,2,3}] [--completeness ] [--host {1,2,3,4}] [--metadata {1,2,3,4}] [--countries ("country1, country2") ] [--remove] [--phylogeny]

Download and filter Lassa virus sequences

options:
  -h, --help            show this help message and exit
  -o                    Output directory for sequences
  --genome {1,2,3}      Genome completeness filter:
                        1 = Complete genomes only (>99 percent of reference length)
                        2 = Partial genomes (specify minimum percent with --completeness)
                        3 = No completeness filter
  --completeness        Minimum sequence completeness (1-100 percent)
                        Required when --genome=2
  --host {1,2,3,4}      Host filter:
                        1 = Human sequences only
                        2 = Rodent sequences only
                        3 = Both human and rodent sequences
                        4 = No host filter
  --metadata {1,2,3,4}  Metadata completeness filter:
                        1 = Keep only sequences with known location
                        2 = Keep only sequences with known date
                        3 = Keep only sequences with both known location and date
                        4 = No metadata filter
  --countries           (Optional) Comma-separated list of countries to filter sequences
                        If not specified, sequences from all countries will be included
                        Examples: "Sierra Leone, Guinea" or "Nigeria, Mali"
                        Available: Nigeria, Sierra Leone, Liberia, Guinea, Mali, Ghana, Benin, Burkina Faso, Ivory Coast, Togo
  --remove              (Optional) File containing accession numbers to remove
                        One accession number per line, lines starting with # are ignored
  --phylogeny          (Optional) Create concatenated FASTA files and perform phylogenetic analysis
                        Includes sequence alignment, trimming, and tree building

example:
  # Download complete genomes from human hosts with known location and date from multiple countries:
  lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --countries "Sierra Leone, Guinea"
```

### Interactive Mode
If optional arguments are not provided, the program will run in interactive mode and prompt for choices:
```bash
lassaseq -o output_directory
```

### Command Line Mode
All options can be specified directly:
```bash
# Download complete genomes from human hosts with known location and date
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3

# Download sequences with ≥80% completeness from both human and rodent hosts
lassaseq -o lassa_output --genome 2 --completeness 80 --host 3 --metadata 4

# Download all sequences without any filtering
lassaseq -o lassa_output --genome 3 --host 4 --metadata 4

# Download sequences from specific countries
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --countries "Sierra Leone, Guinea"

# Download sequences excluding specific accession numbers
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --remove remove.txt

# Download sequences and perform phylogenetic analysis
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --phylogeny
```

### Output Structure

```
output_directory/
├── FASTA/
│   ├── L_segment/
│   │   ├── lassa_l_segments.fasta  (downloaded sequences)
│   │   ├── reference.fasta         (reference sequence NC_004297.1)
│   │   └── outgroup.fasta          (Pinneo strain KM822127.1, Nigeria 1969)
│   ├── S_segment/
│   │   ├── lassa_s_segments.fasta  (downloaded sequences)
│   │   ├── reference.fasta         (reference sequence NC_004296.1)
│   │   └── outgroup.fasta          (Pinneo strain KM822128.1, Nigeria 1969)
│   └── unknown_segment/
│       └── lassa_unknown_segments.fasta
├── Phylogeny/
│   ├── FASTA/
│   │   ├── L_segment/
│   │   │   └── all_l_segments.fasta    (concatenated, deduplicated sequences)
│   │   └── S_segment/
│   │       └── all_s_segments.fasta    (concatenated, deduplicated sequences)
│   ├── MSA/                        (MAFFT alignments)
│   │   ├── L_segment/
│   │   │   └── l_aligned.fasta
│   │   └── S_segment/
│   │       └── s_aligned.fasta
│   ├── TrimAl/                     (Trimmed alignments)
│   │   ├── L_segment/
│   │   │   └── l_trimmed.fasta
│   │   └── S_segment/
│   │       └── s_trimmed.fasta
│   └── Tree/                       (IQ-TREE output)
│       ├── L_segment/
│       │   ├── l_trimmed.fasta.treefile
│       │   ├── l_trimmed.fasta.contree
│       │   └── l_trimmed.fasta.iqtree
│       └── S_segment/
│           ├── s_trimmed.fasta.treefile
│           ├── s_trimmed.fasta.contree
│           └── s_trimmed.fasta.iqtree
└── summary_Lassa.txt
```

### Summary File Content
The summary_Lassa.txt file provides detailed information about:
- Initial dataset statistics
- Detailed host distribution with sequence counts for:
  - Human hosts (e.g., Homo sapiens, human patient)
  - Rodent hosts (e.g., Mastomys natalensis, Hylomyscus pamfi)
  - Other hosts
  - Sequences with no host information
- Geographical distribution at each filtering step
- Impact of completeness filtering (if applied)
- Impact of host filtering (if applied)
- Impact of metadata filtering (if applied)
- Impact of country filtering (if applied)
- Final sequence counts for each segment

### Phylogenetic Analysis
When the `--phylogeny` flag is used, LassaSeq performs the following steps:
1. **Sequence Concatenation**: Combines downloaded sequences with reference and outgroup sequences, removing duplicates.
2. **Multiple Sequence Alignment**: Uses MAFFT with optimized parameters:
   - `--localpair`: For accurate alignment
   - `--maxiterate 1000`: Maximum number of iterative refinement cycles
   - `--ep 0.123`: Gap extension penalty
   - `--op 1.53`: Gap opening penalty
   - `--reorder`: Output order aligned with input order
   - `--adjustdirection`: Adjust sequence directions
3. **Alignment Trimming**: Uses TrimAl with `-automated1` for automated trimming
4. **Tree Building**: Uses IQ-TREE2 with:
   - `-nt AUTO`: Automatic thread optimization
   - `-bb 10000`: Ultra-fast bootstrap with 10,000 replicates
   - `-m TEST`: Automatic model selection

## Requirements
- Python ≥ 3.6
- Biopython
- requests
- MAFFT
- TrimAl
- IQ-TREE2
