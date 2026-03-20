# FunPipe

Automated functional annotation of gene models predicted by [BRAKER](https://github.com/Gaius-Augustus/BRAKER). Starting from a GFF3 and a genome FASTA, the pipeline extracts protein sequences, searches them against **UniProt Swiss-Prot** or **RefSeq**, detects **Pfam** domains, and produces an annotated GenBank file, updated GFF3, and statistics summary.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Databases](#databases)
- [Output Files](#output-files)
- [Pipeline Steps](#pipeline-steps)
- [Citation](#citation)

---

## Overview

```
BRAKER GFF3 + Genome FASTA
         │
         ▼
   [gffread]   Extract & translate CDS → protein FASTA
         │
         ▼
   [blastp]    Search UniProt Swiss-Prot / RefSeq / custom DB
         │
         ▼
   [hmmscan]   Scan against Pfam-A HMMs
         │
         ▼
   [Python]    Parse hits → assign product names & Pfam domains
         │
         ▼
   [BioPython] Build annotated GFF3 + GenBank file
         │
         ▼
  rice_annotation.gbk  +  braker_annotated.gff3  +  annotation_summary.txt
```

The pipeline is split into two scripts that must live in the same directory:

| Script | Role |
|---|---|
| `annotate.sh` | Entry point — argument parsing, dependency checks, database downloads, config generation |
| `functional_annotation_pipeline.sh` | Core engine — runs gffread, BLAST, HMMER, and all Python integration steps |

### How they communicate

`annotate.sh` resolves all input paths to absolute paths, writes them into `annotation_config.sh`, then calls the core script and passes the genome FASTA **explicitly as `$1`**:

```
annotate.sh
  └─ writes annotation_config.sh   (all resolved absolute paths)
  └─ calls functional_annotation_pipeline.sh "$GENOME_FASTA"
                └─ reads $1 → GENOME_FASTA
                └─ sources annotation_config.sh for everything else
```

The genome FASTA is passed as a positional argument rather than being read from the config file. This avoids a shell heredoc expansion issue — variables inside quoted heredocs (`<< 'DELIMITER'`) are not expanded by the shell. All five embedded Python scripts receive their shell-side values as `sys.argv` arguments for the same reason.

---

## Requirements

| Tool | Purpose | Version |
|---|---|---|
| [gffread](https://github.com/gpertea/gffread) | Extract and translate CDS from GFF3 | ≥ 0.12 |
| [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | Protein homology search | ≥ 2.9 |
| [HMMER](http://hmmer.org/) | Pfam domain detection | ≥ 3.3 |
| Python | Result parsing and GenBank generation | ≥ 3.6 |
| [BioPython](https://biopython.org/) | GenBank file construction | ≥ 1.78 |

---

## Installation

**Conda (recommended)**
```bash
conda create -n annotation -c bioconda gffread blast hmmer biopython
conda activate annotation
```

**Ubuntu / Debian**
```bash
sudo apt-get install gffread ncbi-blast+ hmmer python3-biopython
```

**macOS (Homebrew)**
```bash
brew install gffread blast hmmer
pip install biopython
```

**Clone the repository**
```bash
git clone https://github.com/your-username/functional-annotation-pipeline.git
cd functional-annotation-pipeline
chmod +x annotate.sh functional_annotation_pipeline.sh
```

---

## Quick Start

```bash
./annotate.sh -g braker.gff3 -f genome.fasta
```

On first run this will download and format the UniProt Swiss-Prot database (~500 MB) and Pfam-A (~1 GB) automatically before running annotation.

---

## Usage

```
./annotate.sh -g BRAKER.GFF3 -f GENOME.FASTA [OPTIONS]
```

> **Note:** `functional_annotation_pipeline.sh` is called automatically by `annotate.sh` — do not invoke it directly.

### Required arguments

| Flag | Description |
|---|---|
| `-g, --gff FILE` | BRAKER GFF3 file |
| `-f, --fasta FILE` | Genome FASTA file |

### Optional arguments

| Flag | Default | Description |
|---|---|---|
| `-o, --output DIR` | `functional_annotation` | Output directory |
| `-d, --database DB` | `uniprot` | `uniprot` \| `refseq` \| `custom /path/to/proteins.fasta` |
| `-e, --evalue VALUE` | `1e-5` | BLAST e-value threshold |
| `-t, --threads NUM` | `8` | CPU threads for BLAST and HMMER |
| `-c, --config FILE` | auto-generated | Path to a custom configuration file |
| `--setup-only` | — | Download and format databases, then exit |
| `--skip-pfam` | — | Skip Pfam domain detection |
| `--force` | — | Overwrite an existing output directory |
| `-h, --help` | — | Show help |
| `-v, --version` | — | Show version |

### Examples

```bash
# Basic run
./annotate.sh -g braker.gff3 -f genome.fasta

# Custom output directory, 16 threads
./annotate.sh -g braker.gff3 -f genome.fasta -o results/ -t 16

# Use RefSeq plant proteins instead of UniProt
./annotate.sh -g braker.gff3 -f genome.fasta -d refseq

# Use your own protein FASTA as the database
./annotate.sh -g braker.gff3 -f genome.fasta -d /path/to/proteins.fasta

# Download and prepare databases only (no annotation)
./annotate.sh --setup-only

# Skip Pfam domain detection, overwrite previous output
./annotate.sh -g braker.gff3 -f genome.fasta --skip-pfam --force
```

---

## Databases

### UniProt Swiss-Prot (default)
Downloaded automatically on first run.
- Source: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
- Size: ~500 MB compressed
- Manually reviewed, curated protein records with reliable product names

### RefSeq Plant Proteins
Must be downloaded manually (~2–3 GB).
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/plant.protein.fasta.gz
gunzip plant.protein.fasta.gz
mv plant.protein.fasta annotation_setup/databases/refseq_plant.fasta

# Then run with:
./annotate.sh -g braker.gff3 -f genome.fasta -d refseq
```

### Custom Database
Pass any protein FASTA with `-d`. The pipeline formats it with `makeblastdb` automatically.
```bash
./annotate.sh -g braker.gff3 -f genome.fasta -d /path/to/custom_proteins.fasta
```

### Pfam-A
Downloaded automatically on first run.
- Source: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
- Size: ~1 GB compressed
- Indexed with `hmmpress` automatically on first use
- Skip entirely with `--skip-pfam` if domain annotation is not needed

---

## Output Files

All files are written to the output directory (`functional_annotation/` by default).

```
functional_annotation/
├── rice_annotation.gbk       ← Main output (GenBank format)
├── braker_annotated.gff3     ← Input GFF3 enriched with product= and note=Domains: attributes
├── annotation_summary.txt    ← Statistics: annotation rate, top 20 functions, domain coverage
├── blastp_results.txt        ← Raw BLAST tabular output (outfmt 6 + stitle)
├── blastp_stats.txt          ← BLAST summary: mean identity, e-value distribution
├── gene_products.txt         ← Gene ID → product name (TSV)
├── gene_domains.txt          ← Gene ID → Pfam domain(s) (TSV)
├── hmmscan_domains.txt       ← Raw hmmscan domain table (--domtblout)
├── hmmscan_output.txt        ← Full hmmscan stdout
├── braker_proteins.fasta     ← Protein sequences extracted from GFF3 + genome
└── annotation_config.sh      ← Auto-generated config (records all parameters for reproducibility)
```

---

## Pipeline Steps

| Step | Tool | Description |
|---|---|---|
| 1 | `gffread` | Extract CDS from GFF3 + genome FASTA and translate to protein FASTA |
| 2 | `blastp` | Search proteins against chosen database; keep top 5 hits per query (outfmt 6) |
| 3 | `hmmscan` | Scan proteins against Pfam-A HMMs at E ≤ 0.01; skipped with `--skip-pfam` |
| 4 | Python | Parse BLAST tabular output; assign best product name per gene |
| 5 | Python | Parse hmmscan domain table (`--domtblout`); assign Pfam domain list per gene |
| 6 | Python | Write annotated GFF3 with `product=` and `note=Domains:` attributes on CDS features |
| 7 | Python | Generate `annotation_summary.txt` with counts, coverage rates, and top-20 function table |
| 8 | BioPython | Load genome FASTA + annotated GFF3; write final GenBank file |

### Variable passing in embedded Python

All five Python blocks use fully quoted heredocs (`<< 'PYTHON_SCRIPT'`) so the Python source is never touched by the shell. Shell variables are passed as positional arguments and read with `sys.argv`:

```bash
# Shell variables are expanded here, safely outside the heredoc
python3 - "$OUTPUT_DIR/blastp_results.txt" "$OUTPUT_DIR/gene_products.txt" << 'PYTHON_SCRIPT'
import sys
blastp_file, output_file = sys.argv[1], sys.argv[2]
...
PYTHON_SCRIPT
```

---

## Project Structure

```
.
├── annotate.sh                        ← Entry point (run this)
├── functional_annotation_pipeline.sh  ← Core pipeline (called by annotate.sh)
├── annotation_setup/
│   └── databases/
│       ├── uniprot_sprot.fasta
│       ├── uniprot_sprot_db.*         ← BLAST-formatted database files
│       ├── Pfam-A.hmm
│       └── Pfam-A.hmm.h3*            ← HMMER index files
└── functional_annotation/            ← Output directory (created at runtime)
    ├── annotation_config.sh          ← Auto-generated run config
    └── ...                           ← All result files
```

---

## Citation

If you use this pipeline in published work, please cite the underlying tools:

- **gffread**: Pertea G & Pertea M (2020). GFF Utilities: GffRead and GffCompare. *F1000Research* 9:304.
- **BLAST+**: Camacho C et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- **HMMER**: Eddy SR (2011). Accelerated Profile HMM Searches. *PLoS Computational Biology* 7(10):e1002195.
- **UniProt**: UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase. *Nucleic Acids Research* 51:D523–D531.
- **Pfam**: Mistry J et al. (2021). Pfam: The protein families database in 2021. *Nucleic Acids Research* 49:D412–D419.
- **BioPython**: Cock PJA et al. (2009). Biopython: freely available Python tools for computational molecular biology. *Bioinformatics* 25(11):1422–1423.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
