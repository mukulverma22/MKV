#!/bin/bash

################################################################################
# FUNCTIONAL ANNOTATION PIPELINE - User-Friendly Wrapper
# 
# Purpose: Automated functional annotation of BRAKER GFF3 with protein and 
#          domain information from UniProt and Pfam
#
# Usage: ./annotate.sh -g BRAKER.GFF3 -f GENOME.FASTA [OPTIONS]
#
# Author: Your Name
# Version: 1.0.0
# License: MIT
################################################################################

set -e

# ============================================================================
# COLORS FOR OUTPUT
# ============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ============================================================================
# FUNCTIONS
# ============================================================================

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

usage() {
    cat << EOF
${BLUE}FUNCTIONAL ANNOTATION PIPELINE${NC}

${BLUE}USAGE:${NC}
    ./annotate.sh -g BRAKER.GFF3 -f GENOME.FASTA [OPTIONS]

${BLUE}REQUIRED ARGUMENTS:${NC}
    -g, --gff FILE          BRAKER GFF3 file
    -f, --fasta FILE        Genome FASTA file

${BLUE}OPTIONAL ARGUMENTS:${NC}
    -o, --output DIR        Output directory (default: functional_annotation)
    -d, --database DB       Protein database to use:
                            - uniprot (default, recommended)
                            - refseq (NCBI RefSeq plant proteins)
                            - custom PATH (path to FASTA file)
    -e, --evalue VALUE      BLAST e-value threshold (default: 1e-5)
    -t, --threads NUM       Number of CPU threads (default: 8)
    -c, --config FILE       Custom configuration file
    --setup-only            Only setup databases, don't run annotation
    --skip-pfam             Skip Pfam domain detection
    --force                 Overwrite existing output
    -h, --help              Show this help message
    -v, --version           Show version

${BLUE}EXAMPLES:${NC}
    # Basic usage
    ./annotate.sh -g braker.gff3 -f genome.fasta

    # With custom settings
    ./annotate.sh -g braker.gff3 -f genome.fasta -o results -t 16

    # Using RefSeq database
    ./annotate.sh -g braker.gff3 -f genome.fasta -d refseq

    # Setup databases only
    ./annotate.sh --setup-only

${BLUE}OUTPUT FILES:${NC}
    functional_annotation/
    ├── rice_annotation.gbk          ← Main output (GenBank format)
    ├── braker_annotated.gff3        ← Annotated GFF3
    ├── annotation_summary.txt       ← Summary statistics
    ├── blastp_results.txt           ← BLAST hits
    ├── gene_products.txt            ← Product assignments
    ├── gene_domains.txt             ← Pfam domains
    └── braker_proteins.fasta        ← Extracted proteins

${BLUE}REQUIREMENTS:${NC}
    • gffread (part of cufflinks)
    • BLAST+ (blastp, makeblastdb)
    • HMMER (hmmscan, hmmpress)
    • Python 3.6+
    • BioPython

${BLUE}INSTALLATION:${NC}
    # Ubuntu/Debian
    sudo apt-get install gffread ncbi-blast+ hmmer python3-biopython

    # macOS with Homebrew
    brew install gffread blast hmmer
    pip install biopython

    # Conda
    conda create -n annotation -c bioconda gffread blast hmmer biopython
    conda activate annotation

EOF
    exit 0
}

version() {
    echo "Functional Annotation Pipeline v1.0.0"
    exit 0
}

check_requirements() {
    print_info "Checking requirements..."
    
    local missing=0
    
    local commands=("gffread" "blastp" "makeblastdb" "hmmscan" "hmmpress" "python3")
    
    for cmd in "${commands[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            print_error "$cmd not found"
            missing=$((missing + 1))
        fi
    done
    
    if [ $missing -gt 0 ]; then
        echo ""
        print_error "Missing $missing required tool(s)"
        echo ""
        echo "Install with:"
        echo "  Ubuntu/Debian: sudo apt-get install gffread ncbi-blast+ hmmer"
        echo "  macOS: brew install gffread blast hmmer"
        echo "  Conda: conda install -c bioconda gffread blast hmmer"
        exit 1
    fi
    
    python3 -c "import Bio" 2>/dev/null || {
        print_error "BioPython not found"
        echo "Install with: pip install biopython"
        exit 1
    }
    
    print_success "All requirements satisfied"
}

check_input_files() {
    if [ ! -f "$GFF_FILE" ]; then
        print_error "GFF3 file not found: $GFF_FILE"
        exit 1
    fi
    
    if [ ! -f "$FASTA_FILE" ]; then
        print_error "FASTA file not found: $FASTA_FILE"
        exit 1
    fi
    
    print_success "Input files verified"
}

setup_databases() {
    print_header "Setting Up Databases"
    
    mkdir -p "$DB_DIR"
    
    case "$DATABASE" in
        uniprot)
            setup_uniprot
            PROTEIN_DB="$DB_DIR/uniprot_sprot_db"
            ;;
        refseq)
            setup_refseq
            PROTEIN_DB="$DB_DIR/refseq_plant_db"
            ;;
        *)
            # Custom database path
            if [ ! -f "$DATABASE" ]; then
                print_error "Custom database file not found: $DATABASE"
                exit 1
            fi
            CUSTOM_DB=$(basename "$DATABASE" .fasta)
            cp "$DATABASE" "$DB_DIR/$CUSTOM_DB.fasta"
            makeblastdb -in "$DB_DIR/$CUSTOM_DB.fasta" -dbtype prot \
                -out "$DB_DIR/${CUSTOM_DB}_db" -title "$CUSTOM_DB" 2>/dev/null
            PROTEIN_DB="$DB_DIR/${CUSTOM_DB}_db"
            ;;
    esac
    
    if [ "$SKIP_PFAM" -eq 0 ]; then
        setup_pfam
    fi
    
    print_success "Database setup complete"
}

setup_uniprot() {
    if [ -f "$DB_DIR/uniprot_sprot.fasta" ] && [ -s "$DB_DIR/uniprot_sprot.fasta" ]; then
        print_info "UniProt database already exists"
    else
        print_info "Downloading UniProt Swiss-Prot..."
        print_warning "This may take 5-10 minutes (~500 MB)"
        
        cd "$DB_DIR"
        
        if curl -L --max-time 3600 \
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
            -o uniprot_sprot.fasta.gz 2>/dev/null; then
            
            gunzip -f uniprot_sprot.fasta.gz
            print_success "UniProt downloaded"
        else
            print_error "Failed to download UniProt"
            print_info "Manual download required:"
            echo "  1. Visit: https://www.uniprot.org/downloads"
            echo "  2. Download: uniprot_sprot.fasta.gz"
            echo "  3. Extract: gunzip uniprot_sprot.fasta.gz"
            echo "  4. Copy to: $DB_DIR/uniprot_sprot.fasta"
            exit 1
        fi
        
        cd - > /dev/null
    fi
    
    if [ ! -f "$DB_DIR/uniprot_sprot_db.phr" ]; then
        print_info "Formatting UniProt for BLAST..."
        makeblastdb -in "$DB_DIR/uniprot_sprot.fasta" -dbtype prot \
            -out "$DB_DIR/uniprot_sprot_db" -title "UniProt SwissProt" 2>/dev/null
    fi
}

setup_refseq() {
    if [ -f "$DB_DIR/refseq_plant.fasta" ] && [ -s "$DB_DIR/refseq_plant.fasta" ]; then
        print_info "RefSeq database already exists"
    else
        print_error "RefSeq database not found"
        echo ""
        echo "RefSeq plant proteins must be downloaded manually (2-3 GB)"
        echo "Visit: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/"
        echo "Download: plant.protein.fasta.gz"
        echo ""
        echo "Then:"
        echo "  gunzip plant.protein.fasta.gz"
        echo "  mv plant.protein.fasta $DB_DIR/refseq_plant.fasta"
        exit 1
    fi
    
    if [ ! -f "$DB_DIR/refseq_plant_db.phr" ]; then
        print_info "Formatting RefSeq for BLAST..."
        makeblastdb -in "$DB_DIR/refseq_plant.fasta" -dbtype prot \
            -out "$DB_DIR/refseq_plant_db" -title "RefSeq Plant" 2>/dev/null
    fi
}

setup_pfam() {
    if [ -f "$DB_DIR/Pfam-A.hmm" ] && [ -s "$DB_DIR/Pfam-A.hmm" ]; then
        if [ -f "$DB_DIR/Pfam-A.hmm.h3m" ]; then
            print_info "Pfam database already exists and is indexed"
            return
        fi
    else
        print_info "Downloading Pfam-A..."
        print_warning "This may take 10-15 minutes (~1 GB)"
        
        cd "$DB_DIR"
        
        if curl -L --max-time 3600 \
            "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" \
            -o Pfam-A.hmm.gz 2>/dev/null; then
            
            gunzip -f Pfam-A.hmm.gz
            print_success "Pfam downloaded"
        else
            print_warning "Pfam download failed, will skip domain detection"
            SKIP_PFAM=1
            cd - > /dev/null
            return
        fi
        
        cd - > /dev/null
    fi
    
    if [ ! -f "$DB_DIR/Pfam-A.hmm.h3m" ]; then
        print_info "Indexing Pfam database (this may take a few minutes)..."
        hmmpress -f "$DB_DIR/Pfam-A.hmm" 2>/dev/null
    fi
}

create_config_file() {
    # Resolve absolute paths before writing to config so the core
    # pipeline script can be run from any working directory.
    local abs_gff abs_fasta
    abs_gff="$(cd "$(dirname "$GFF_FILE")" && pwd)/$(basename "$GFF_FILE")"
    abs_fasta="$(cd "$(dirname "$FASTA_FILE")" && pwd)/$(basename "$FASTA_FILE")"

    cat > "$CONFIG_FILE" << EOF
#!/bin/bash

# Auto-generated configuration file
# Generated: $(date)

export BRAKER_GFF3="${abs_gff}"
export GENOME_FASTA="${abs_fasta}"
export PROTEIN_DB="${PROTEIN_DB}"
export PROTEIN_FASTA="${PROTEIN_FASTA}"
export PFAM_DB="${DB_DIR}/Pfam-A.hmm"
export BLAST_EVALUE="${EVALUE}"
export BLAST_MAX_TARGETS="5"
export THREADS="${THREADS}"
export OUTPUT_DIR="${OUTPUT_DIR}"
export GFFREAD_BIN="gffread"
export BLASTP_BIN="blastp"
export MAKEBLASTDB_BIN="makeblastdb"
export HMMSCAN_BIN="hmmscan"
export HMMPRESS_BIN="hmmpress"
export SKIP_PFAM="${SKIP_PFAM}"

EOF
}

run_annotation() {
    print_header "Running Functional Annotation"
    
    mkdir -p "$OUTPUT_DIR"
    
    # Source config and run the core pipeline, passing GENOME_FASTA
    # explicitly as a command-line argument to avoid heredoc expansion issues.
    source "$CONFIG_FILE"
    
    print_info "GFF3 Input: $BRAKER_GFF3"
    print_info "Genome:     $GENOME_FASTA"
    print_info "Protein DB: $PROTEIN_DB"
    print_info "Threads:    $THREADS"
    print_info "Output:     $OUTPUT_DIR"
    
    bash "${SCRIPT_DIR}/functional_annotation_pipeline.sh" "$GENOME_FASTA" 2>&1
    
    print_success "Annotation complete!"
}

show_summary() {
    echo ""
    print_header "Summary"
    
    if [ -f "$OUTPUT_DIR/annotation_summary.txt" ]; then
        cat "$OUTPUT_DIR/annotation_summary.txt"
    fi
    
    echo ""
    print_success "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Main output file:"
    echo "  → $OUTPUT_DIR/rice_annotation.gbk"
    echo ""
}

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Default values
GFF_FILE=""
FASTA_FILE=""
OUTPUT_DIR="functional_annotation"
DATABASE="uniprot"
EVALUE="1e-5"
THREADS="8"
SETUP_ONLY=0
SKIP_PFAM=0
FORCE=0
CONFIG_FILE=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--gff)
            GFF_FILE="$2"
            shift 2
            ;;
        -f|--fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--database)
            DATABASE="$2"
            shift 2
            ;;
        -e|--evalue)
            EVALUE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --setup-only)
            SETUP_ONLY=1
            shift
            ;;
        --skip-pfam)
            SKIP_PFAM=1
            shift
            ;;
        --force)
            FORCE=1
            rm -rf "$OUTPUT_DIR"
            shift
            ;;
        -h|--help)
            usage
            ;;
        -v|--version)
            version
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$GFF_FILE" ] || [ -z "$FASTA_FILE" ]; then
    if [ "$SETUP_ONLY" -eq 0 ]; then
        print_error "Missing required arguments: -g (GFF3) and -f (FASTA) are required"
        echo ""
        usage
    fi
fi

# Set default config file path
if [ -z "$CONFIG_FILE" ]; then
    CONFIG_FILE="$OUTPUT_DIR/annotation_config.sh"
fi

# Main execution
print_header "FUNCTIONAL ANNOTATION PIPELINE v1.0.0"

# Check if output already exists
if [ -d "$OUTPUT_DIR" ] && [ "$FORCE" -eq 0 ]; then
    print_warning "Output directory already exists: $OUTPUT_DIR"
    read -p "Overwrite? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Exiting"
        exit 0
    fi
    rm -rf "$OUTPUT_DIR"
fi

# Check requirements
check_requirements

# Check input files (skip if setup-only)
if [ "$SETUP_ONLY" -eq 0 ]; then
    check_input_files
fi

# Create database directory
DB_DIR="$(pwd)/annotation_setup/databases"

# Setup databases
setup_databases

# Create config file
mkdir -p "$OUTPUT_DIR"
create_config_file

print_success "Configuration saved to: $CONFIG_FILE"

# Run annotation if not setup-only
if [ "$SETUP_ONLY" -eq 0 ]; then
    echo ""
    run_annotation
    show_summary
else
    print_success "Database setup complete!"
    print_info "To run annotation:"
    echo "  ./annotate.sh -g $GFF_FILE -f $FASTA_FILE"
fi

exit 0
