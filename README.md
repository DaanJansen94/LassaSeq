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

**Important**: Always activate the conda environment before running any LassaSeq commands:
```bash
conda activate lassa_env
```

To see all available options:
```bash
lassaseq --help
```

This will display:
```
usage: lassaseq [-h] -o  [--genome {1,2,3}] [--completeness ] [--host {1,2,3,4}] [--metadata {1,2,3,4}] [--countries ("country1, country2") ] [--remove] [--phylogeny] [--consensus_L] [--consensus_S] [--lineage LINEAGE] [--sublineage SUBLINEAGE]

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
  --phylogeny           (Optional) Create concatenated FASTA files and perform phylogenetic analysis
                        Includes sequence alignment, trimming, and tree building
  --consensus_L         (Optional) Path to custom consensus sequences for L segment
  --consensus_S         (Optional) Path to custom consensus sequences for S segment
  --lineage             (Optional) Filter sequences by lineage
  --sublineage          (Optional) Filter sequences by sublineage                    
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

# Download sequences from specific countries and lineage
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --countries "Sierra Leone, Guinea" --lineage IV

# Download sequences from specific lineage and sublineage
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --phylogeny --lineage IV --sublineage III

# Download sequences and perform phylogenetic analysis with lineage filtering
lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --phylogeny --lineage IV
```
### Output Structure

```
output_directory/
├── FASTA/
│   ├── L_segment/
│   │   ├── lassa_l_segments.fasta  (downloaded sequences)
│   │   ├── reference.fasta         (reference sequence NC_004297.1)
│   │   ├── outgroup.fasta          (Pinneo strain KM822127.1, Nigeria 1969)
│   │   └── consensus.fasta         (optional custom consensus sequence(s))
│   ├── S_segment/
│   │   ├── lassa_s_segments.fasta  (downloaded sequences)
│   │   ├── reference.fasta         (reference sequence NC_004296.1)
│   │   ├── outgroup.fasta          (Pinneo strain KM822128.1, Nigeria 1969)
│   │   └── consensus.fasta         (optional custom consensus sequence(s))
│   └── unknown_segment/
│       └── lassa_unknown_segments.fasta
├── Phylogeny/
│   ├── FASTA/
│   │   ├── L_segment/
│   │   │   ├── all_l_segments.fasta    (concatenated, deduplicated sequences)
│   │   │   └── l_metadata.txt          (metadata for FigTree visualization)
│   │   └── S_segment/
│   │       ├── all_s_segments.fasta    (concatenated, deduplicated sequences)
│   │       └── s_metadata.txt          (metadata for FigTree visualization)
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

### Lineage Information
The tool includes lineage information for Lassa virus sequences in the `lineages` directory. The lineages are in line with current literature, and sublineages for lineage IV were determined using fastbaps software. Users can modify or update the lineage files (`l_lineages.txt` and `s_lineages.txt`) to include their own lineage assignments. Each file follows a tab-delimited format with columns for accession numbers, lineages, and sublineages.

### Summary File Content
The summary_Lassa.txt file provides detailed information about:
- Detailed host distribution with sequence counts for:
  - Human hosts (e.g., Homo sapiens, human patient)
  - Rodent hosts (e.g., Mastomys natalensis, Hylomyscus pamfi)
  - Other hosts
  - Sequences with no host information
- Geographical distribution at each filtering step
- Final sequence counts for each segment

### Phylogenetic Analysis
When the `--phylogeny` flag is used, LassaSeq performs the following steps:
1. **Sequence Concatenation**: Combines downloaded sequences with reference and outgroup sequences, removing duplicates.
2. **Multiple Sequence Alignment**: Uses MAFFT with optimized parameters:
3. **Alignment Trimming**: Uses TrimAl with `-automated1` for automated trimming
4. **Tree Building**: Uses IQ-TREE2 with ultra-fast bootstrap and 10,000 replicates

### Tree Visualization
The phylogenetic trees generated by LassaSeq (`.treefile` files) can be visualuzed using various software packages (e.g., FigTree) and annotated using the metadata files provided by LassaSeq.

1. **Install FigTree**:
   - Download from: http://tree.bio.ed.ac.uk/software/figtree/
   - Available for Windows, Mac, and Linux

2. **Visualizing Trees**:
   - Open FigTree
   - Load the tree file (e.g., `l_trimmed.fasta.treefile` or `s_trimmed.fasta.treefile`)
   - Import the corresponding metadata file (`l_metadata.txt` or `s_metadata.txt`)
   - Use the metadata columns to color and annotate your tree:
     - Location: Country of origin (e.g., SierraLeone, Nigeria)
     - Location2: City or specific location (e.g., Kenema, Nzerekore, Unknown)
     - Host: Host species (Human, Rodent, Other)
     - Date: Collection date in decimal format (e.g., 1969.000)

## Handling Recombination Events
Due to the computationally intensive nature of recombination detection, this analysis is not included in LassaSeq itself. However, it is highly recommended to perform recombination analysis as a separate step before conducting phylogenetic analyses, as recombinant sequences can interfere with accurate tree reconstruction. Here's the recommended workflow:

1. **Recombination Analysis**:
   Before running LassaSeq's phylogenetic analysis, perform recombination detection using your preferred tool. For example, HYPHY GARD can be installed and used as follows:
   ```bash
   # Example using HYPHY GARD
   conda install bioconda::hyphy
   hyphy gard --rv Gamma --mode faster --input trimmed_alignment.fasta --output gard_output.json
   ```
   Note: Recombination analysis can be time-consuming depending on your dataset size and chosen tool.
   
2. **Remove Recombinant Sequences**:
   Once you've identified recombinant sequences using your chosen analysis tool, add their accession numbers to a `remove.txt` file. Use this file with LassaSeq's `--remove` option to exclude these sequences from the analysis:
   ```bash
   lassaseq -o lassa_output --genome 1 --host 1 --metadata 3 --remove remove.txt --phylogeny
   ```

This optional but recommended step is crucial for obtaining reliable phylogenetic trees, as recombinant sequences can lead to incorrect evolutionary relationships and branch patterns.

## Dependencies

- Python ≥ 3.6
- BioPython ≥ 1.79
- NumPy ≥ 1.21.0
- MAFFT
- TrimAl
- IQTree2

## Citation

If you use Lassaseq in your research, please cite:

```
Jansen, D., & Vercauteren, K. (2025). LassaSeq: A Command-Line Tool for Downloading, Processing and Analyzing Lassa Virus Sequences for Phylogenetic Analysis (Version v0.1.1) [Computer software]. https://doi.org/10.5281/zenodo.14936276
```

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.
