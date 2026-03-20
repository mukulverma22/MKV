#!/bin/bash

# FunPipe for GFF3 
# Updated for: UniProt Swiss-Prot + Pfam-A databases
#
# Usage (called by annotate.sh):
#   functional_annotation_pipeline.sh <GENOME_FASTA>
#
# GENOME_FASTA is passed as $1 to avoid shell variable expansion issues
# inside Python heredocs. All other parameters come from annotation_config.sh.

set -e

# ============================================================================
# ARGUMENT + CONFIG
# ============================================================================

# $1 — absolute path to genome FASTA, passed explicitly by annotate.sh
GENOME_FASTA_ARG="$1"

if [ -z "$GENOME_FASTA_ARG" ]; then
    echo "Error: genome FASTA path must be passed as the first argument."
    echo "Usage: functional_annotation_pipeline.sh <GENOME_FASTA>"
    exit 1
fi

if [ ! -f "$GENOME_FASTA_ARG" ]; then
    echo "Error: Genome FASTA file not found: $GENOME_FASTA_ARG"
    exit 1
fi

# Load remaining configuration (set by annotate.sh)
CONFIG_FILE="annotation_config.sh"
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "Error: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Override GENOME_FASTA with the explicitly passed argument so every step
# uses the same resolved path regardless of what the config contains.
GENOME_FASTA="$GENOME_FASTA_ARG"

# ============================================================================
# VALIDATE INPUTS
# ============================================================================

if [ ! -f "$BRAKER_GFF3" ]; then
    echo "Error: BRAKER GFF3 file not found: $BRAKER_GFF3"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "FUNCTIONAL ANNOTATION PIPELINE FOR RICE GENOME"
echo "=========================================="
echo "Input BRAKER GFF3: $BRAKER_GFF3"
echo "Genome FASTA:      $GENOME_FASTA"
echo "Protein DB:        $PROTEIN_DB"
echo "Output directory:  $OUTPUT_DIR"
echo "=========================================="
echo ""

# ============================================================================
# STEP 1: Extract predicted proteins from BRAKER GFF3
# ============================================================================
echo "[Step 1] Extracting predicted proteins from BRAKER GFF3..."

$GFFREAD_BIN "$BRAKER_GFF3" -g "$GENOME_FASTA" -y "$OUTPUT_DIR/braker_proteins.fasta"

PROTEIN_COUNT=$(grep -c "^>" "$OUTPUT_DIR/braker_proteins.fasta")
echo "✓ Proteins extracted: $PROTEIN_COUNT sequences"
echo "  File: $OUTPUT_DIR/braker_proteins.fasta"
echo ""

# ============================================================================
# STEP 2: BLASTP against UniProt Swiss-Prot (or chosen database)
# ============================================================================
echo "[Step 2] Running BLASTP against protein database..."
echo "  Parameters: evalue=$BLAST_EVALUE, max_targets=$BLAST_MAX_TARGETS, threads=$THREADS"

$BLASTP_BIN \
    -query  "$OUTPUT_DIR/braker_proteins.fasta" \
    -db     "$PROTEIN_DB" \
    -out    "$OUTPUT_DIR/blastp_results.txt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -evalue "$BLAST_EVALUE" \
    -num_threads "$THREADS" \
    -max_target_seqs "$BLAST_MAX_TARGETS"

BLASTP_HITS=$(wc -l < "$OUTPUT_DIR/blastp_results.txt")
echo "✓ BLASTP completed with $BLASTP_HITS hits"
echo "  File: $OUTPUT_DIR/blastp_results.txt"
echo ""

# ============================================================================
# STEP 3: hmmscan for Pfam domain identification
# ============================================================================
if [ "$SKIP_PFAM" -eq 0 ]; then
    echo "[Step 3] Running hmmscan for Pfam domain detection..."
    echo "  Database: $PFAM_DB"

    $HMMSCAN_BIN \
        --domtblout "$OUTPUT_DIR/hmmscan_domains.txt" \
        -E 0.01 \
        --cpu "$THREADS" \
        "$PFAM_DB" \
        "$OUTPUT_DIR/braker_proteins.fasta" > "$OUTPUT_DIR/hmmscan_output.txt" 2>&1

    PFAM_HITS=$(grep -v "^#" "$OUTPUT_DIR/hmmscan_domains.txt" | wc -l)
    echo "✓ hmmscan completed with $PFAM_HITS Pfam domain hits"
    echo "  Files: hmmscan_domains.txt, hmmscan_output.txt"
else
    echo "[Step 3] Skipping Pfam domain detection (--skip-pfam)"
    # Write empty placeholder files so downstream steps do not break
    touch "$OUTPUT_DIR/hmmscan_domains.txt"
    touch "$OUTPUT_DIR/hmmscan_output.txt"
fi
echo ""

# ============================================================================
# STEP 4: Parse BLASTP results — extract best product name per gene
# ============================================================================
echo "[Step 4] Parsing BLASTP results for product assignment..."

# Pass shell variables as script arguments; heredoc stays fully quoted so
# no accidental expansion occurs inside the Python code.
python3 - \
    "$OUTPUT_DIR/blastp_results.txt" \
    "$OUTPUT_DIR/gene_products.txt" \
    "$OUTPUT_DIR/blastp_stats.txt" \
<< 'PYTHON_SCRIPT'
import sys
import re

blastp_file, output_file, stats_file = sys.argv[1], sys.argv[2], sys.argv[3]

products      = {}
evalue_list   = []
identity_list = []
coverage_list = {}

with open(blastp_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 13:
            continue

        query_id      = fields[0]
        subject_title = fields[12]
        evalue        = float(fields[10])
        pident        = float(fields[2])
        length        = int(fields[3])

        # Keep only the best (first) hit per query
        if query_id not in products:
            # UniProt title format: sp|ACCESSION|NAME Description [Organism]
            # Take everything before the first '[' and strip the accession token
            clean = re.sub(r'\[.*?\]', '', subject_title).strip()
            parts = clean.split()
            # Drop the leading sp|...|... token if present
            if parts and '|' in parts[0]:
                parts = parts[1:]
            protein_name = ' '.join(parts).strip()
            protein_name = re.sub(r'\s+', ' ', protein_name)

            # GenBank /product qualifier limit
            if len(protein_name) > 255:
                protein_name = protein_name[:255]

            products[query_id] = protein_name
            evalue_list.append(evalue)
            identity_list.append(pident)
            coverage_list[query_id] = length

with open(output_file, 'w') as f:
    for gene_id, product in sorted(products.items()):
        f.write(f"{gene_id}\t{product}\n")

with open(stats_file, 'w') as f:
    n = len(products)
    f.write("BLASTP ANNOTATION STATISTICS\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Total genes with BLAST hits: {n}\n")
    if n:
        f.write(f"Average identity:  {sum(identity_list)/n:.2f}%\n")
        f.write(f"Average coverage:  {sum(coverage_list.values())/n:.0f} aa\n")
        f.write(f"Best e-value:      {min(evalue_list):.2e}\n")
        f.write(f"Worst e-value:     {max(evalue_list):.2e}\n")
        f.write(f"Mean e-value:      {sum(evalue_list)/n:.2e}\n")

print(f"Parsed {len(products)} gene products from BLASTP")
PYTHON_SCRIPT

PRODUCTS_COUNT=$(wc -l < "$OUTPUT_DIR/gene_products.txt")
echo "✓ Product names extracted: $PRODUCTS_COUNT genes"
echo "  File: $OUTPUT_DIR/gene_products.txt"
echo ""

# ============================================================================
# STEP 5: Parse hmmscan results — extract Pfam domains per gene
# ============================================================================
echo "[Step 5] Extracting Pfam domain annotations..."

python3 - \
    "$OUTPUT_DIR/hmmscan_domains.txt" \
    "$OUTPUT_DIR/gene_domains.txt" \
<< 'PYTHON_SCRIPT'
import sys

hmmscan_file, output_file = sys.argv[1], sys.argv[2]

domains = {}

with open(hmmscan_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue

        fields = line.split()
        if len(fields) < 7:
            continue

        # domtblout columns (0-based):
        #  0  target name (query protein)
        #  1  target accession
        #  2  query name  (domain)  — NOTE: hmmscan swaps query/target vs hmmsearch
        #  4  accession of domain
        #  6  sequence E-value
        query_id         = fields[0]
        domain_name      = fields[2]
        domain_accession = fields[4]
        evalue           = float(fields[6])

        if query_id not in domains:
            domains[query_id] = []

        domains[query_id].append(f"{domain_name} ({domain_accession}, e={evalue:.2e})")

with open(output_file, 'w') as f:
    for gene_id in sorted(domains.keys()):
        f.write(f"{gene_id}\t{'; '.join(domains[gene_id])}\n")

print(f"Extracted domains for {len(domains)} genes")
PYTHON_SCRIPT

DOMAIN_COUNT=$(wc -l < "$OUTPUT_DIR/gene_domains.txt")
echo "✓ Pfam domains extracted: $DOMAIN_COUNT genes"
echo "  File: $OUTPUT_DIR/gene_domains.txt"
echo ""

# ============================================================================
# STEP 6: Update GFF3 with product and domain information
# ============================================================================
echo "[Step 6] Updating GFF3 with functional annotations..."

python3 - \
    "$BRAKER_GFF3" \
    "$OUTPUT_DIR/gene_products.txt" \
    "$OUTPUT_DIR/gene_domains.txt" \
    "$OUTPUT_DIR/braker_annotated.gff3" \
<< 'PYTHON_SCRIPT'
import sys
import re

braker_gff, products_file, domains_file, output_gff = \
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

# Load products
products = {}
with open(products_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2:
            products[parts[0]] = parts[1]

# Load domains
domains = {}
with open(domains_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2:
            domains[parts[0]] = parts[1]

cds_count      = 0
annotated_count = 0

with open(braker_gff, 'r') as f_in, open(output_gff, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'):
            f_out.write(line)
            continue

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9 or fields[2] != 'CDS':
            f_out.write(line)
            continue

        cds_count  += 1
        attributes  = fields[8]

        parent_match = re.search(r'Parent=([^;]+)', attributes)
        if parent_match:
            parent_id = parent_match.group(1)

            # Add / replace product qualifier
            if parent_id in products:
                product = products[parent_id]
                if 'product=' not in attributes:
                    attributes += f';product={product}'
                else:
                    attributes = re.sub(r'product=[^;]*', f'product={product}', attributes)
                annotated_count += 1

            # Add domain note (GFF3-escape semicolons and equals signs)
            if parent_id in domains:
                domain_str = domains[parent_id] \
                    .replace(';', '%3B') \
                    .replace('=', '%3D')
                if 'note=' not in attributes:
                    attributes += f';note=Domains: {domain_str}'

        fields[8] = attributes
        f_out.write('\t'.join(fields) + '\n')

print(f"✓ Processed {cds_count} CDS features")
print(f"✓ Annotated {annotated_count} genes with products")
PYTHON_SCRIPT

echo "✓ Annotated GFF3 created: $OUTPUT_DIR/braker_annotated.gff3"
echo ""

# ============================================================================
# STEP 7: Generate summary report
# ============================================================================
echo "[Step 7] Generating summary report..."

python3 - \
    "$OUTPUT_DIR/braker_annotated.gff3" \
    "$OUTPUT_DIR/annotation_summary.txt" \
<< 'PYTHON_SCRIPT'
import sys
import re

gff_file, report_file = sys.argv[1], sys.argv[2]

genes_with_product = 0
genes_with_domains = 0
total_genes        = 0
products_assigned  = {}

with open(gff_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != 'CDS':
            continue

        total_genes += 1
        attributes   = fields[8]

        product_match = re.search(r'product=([^;]+)', attributes)
        if product_match:
            genes_with_product += 1
            product = product_match.group(1)
            products_assigned[product] = products_assigned.get(product, 0) + 1

        if 'note=Domains:' in attributes:
            genes_with_domains += 1

with open(report_file, 'w') as f:
    f.write("=" * 60 + "\n")
    f.write("FUNCTIONAL ANNOTATION SUMMARY REPORT\n")
    f.write("Rice Genome - BRAKER + UniProt + Pfam\n")
    f.write("=" * 60 + "\n\n")

    f.write("OVERVIEW\n")
    f.write("-" * 60 + "\n")
    f.write(f"Total CDS features:                     {total_genes}\n")
    f.write(f"Genes with functional annotation (BLASTP): {genes_with_product}\n")
    if total_genes:
        f.write(f"Annotation coverage:                    {(genes_with_product/total_genes)*100:.2f}%\n")
    f.write(f"Genes with Pfam domains:                {genes_with_domains}\n")
    if total_genes:
        f.write(f"Domain coverage:                        {(genes_with_domains/total_genes)*100:.2f}%\n")

    f.write("\nTOP 20 ASSIGNED FUNCTIONS\n")
    f.write("-" * 60 + "\n")
    for i, (product, count) in enumerate(
            sorted(products_assigned.items(), key=lambda x: x[1], reverse=True)[:20], 1):
        f.write(f"{i:2d}. {product}: {count} genes\n")

    f.write("\nOUTPUT FILES GENERATED\n")
    f.write("-" * 60 + "\n")
    f.write("braker_annotated.gff3  - GFF3 with functional annotations\n")
    f.write("braker_proteins.fasta  - Predicted proteins from BRAKER\n")
    f.write("blastp_results.txt     - BLASTP alignment results\n")
    f.write("blastp_stats.txt       - BLASTP statistics\n")
    f.write("gene_products.txt      - Extracted product names\n")
    f.write("gene_domains.txt       - Pfam domain annotations\n")
    f.write("hmmscan_domains.txt    - Raw hmmscan output\n")
    f.write("rice_annotation.gbk    - GenBank format (final output)\n")

print("Summary report generated")
PYTHON_SCRIPT

cat "$OUTPUT_DIR/annotation_summary.txt"
echo ""

# ============================================================================
# STEP 8: Convert annotated GFF3 to GenBank format
# ============================================================================
echo "[Step 8] Converting annotated GFF3 to GenBank format..."

# GENOME_FASTA is passed as $1 (now stored in GENOME_FASTA_ARG / GENOME_FASTA).
# We pass it explicitly as a script argument so the Python heredoc stays
# fully quoted — no variable expansion risk.
python3 - \
    "$GENOME_FASTA" \
    "$OUTPUT_DIR/braker_annotated.gff3" \
    "$OUTPUT_DIR/rice_annotation.gbk" \
<< 'PYTHON_SCRIPT'
import sys
from datetime import datetime
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

genome_file, gff_file, output_gbk = sys.argv[1], sys.argv[2], sys.argv[3]

print("Reading genome sequences...")
sequences = {}
for record in SeqIO.parse(genome_file, "fasta"):
    sequences[record.id] = record

print(f"Loaded {len(sequences)} sequences from genome")

print("Parsing GFF3 and creating features...")
gene_features = {}

with open(gff_file, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            continue

        seqid        = fields[0]
        feature_type = fields[2]
        start        = int(fields[3]) - 1   # GFF3 is 1-based → 0-based
        end          = int(fields[4])
        strand       = fields[6]
        attributes   = fields[8]

        if seqid not in gene_features:
            gene_features[seqid] = []

        # Parse attributes into a dict; unescape GFF3 percent-encoding
        attr_dict = {}
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                value = value.replace('%3B', ';').replace('%3D', '=')
                attr_dict[key] = value

        strand_val = 1 if strand == '+' else -1
        location   = FeatureLocation(start, end, strand=strand_val)

        qualifiers = {}

        locus_tag = attr_dict.get('locus_tag', attr_dict.get('ID', ''))
        if locus_tag:
            qualifiers['locus_tag'] = [locus_tag]

        qualifiers['product'] = [attr_dict.get('product', 'hypothetical protein')]

        if 'note' in attr_dict:
            qualifiers['note'] = [attr_dict['note']]

        gene_features[seqid].append(
            SeqFeature(location, type=feature_type, qualifiers=qualifiers)
        )

print("Writing GenBank file...")
records = []
for seqid in sorted(sequences.keys()):
    record = sequences[seqid]
    if seqid in gene_features:
        record.features = sorted(gene_features[seqid],
                                 key=lambda feat: feat.location.start)

    record.annotations['organism']      = 'Oryza sativa'
    record.annotations['taxonomy']      = [
        'Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta',
        'Tracheophyta', 'Spermatophyta', 'Magnoliopsida', 'Monocots',
        'Poaceae', 'BEP clade', 'Oryzoideae', 'Oryza'
    ]
    record.annotations['molecule_type'] = 'DNA'
    record.annotations['topology']      = 'linear'
    record.annotations.setdefault(
        'date', datetime.now().strftime('%d-%b-%Y').upper()
    )
    records.append(record)

SeqIO.write(records, output_gbk, "genbank")

total_features = sum(len(r.features) for r in records)
print(f"✓ GenBank file generated: {output_gbk}")
print(f"✓ Total sequences:         {len(records)}")
print(f"✓ Total annotated features:{total_features}")
PYTHON_SCRIPT

echo "✓ GenBank file successfully created"
echo ""

# ============================================================================
# FINAL SUMMARY
# ============================================================================
echo "=========================================="
echo "ANNOTATION PIPELINE COMPLETED!"
echo "=========================================="
echo ""
echo "KEY OUTPUT FILES:"
echo "  📁 $OUTPUT_DIR/"
echo "  ├── rice_annotation.gbk        (✓ MAIN OUTPUT - GenBank format)"
echo "  ├── braker_annotated.gff3      (GFF3 with annotations)"
echo "  ├── annotation_summary.txt     (Summary statistics)"
echo "  ├── blastp_results.txt         (BLASTP hits)"
echo "  ├── blastp_stats.txt           (BLASTP statistics)"
echo "  ├── gene_products.txt          (Product assignments)"
echo "  ├── gene_domains.txt           (Pfam domains)"
echo "  └── hmmscan_output.txt         (hmmscan results)"
echo ""
echo "ANNOTATION STATISTICS:"
echo "  • Protein homology search: UniProt Swiss-Prot"
echo "  • Domain detection: Pfam-A"
echo "  • E-value threshold: $BLAST_EVALUE"
echo ""
echo "✓ Ready for downstream analysis or NCBI submission!"
echo "=========================================="
