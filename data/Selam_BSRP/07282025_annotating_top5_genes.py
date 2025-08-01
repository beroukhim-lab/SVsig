import sys

if len(sys.argv) != 5:
    print("Usage: python 07282025_PCAWG_top5_genes.py sv_file.txt refgene.txt cosmic.csv output_file.tsv")
    sys.exit(1)

sv_path = sys.argv[1]
refgene_path = sys.argv[2]
cosmic_path = sys.argv[3]
output_path = sys.argv[4]

# === Load SV file ===
with open(sv_path, 'r') as f:
    lines = [line.strip().split('\t') for line in f if line.strip()]
    sig_driv_header = lines[0]
    sig_driv_rows = lines[1:]

# === Load refgene ===
with open(refgene_path, 'r') as f:
    refgene_lines = [line.strip().split('\t') for line in f if line.strip()]
    refgene_header = refgene_lines[0]
    refgene_data = refgene_lines[1:]

chr_idx = refgene_header.index('chr')
start_idx = refgene_header.index('start')
end_idx = refgene_header.index('end')
strand_idx = refgene_header.index('strand')
symb_idx = refgene_header.index('symb')
gene_idx = refgene_header.index('gene')

refgene_normalized = []
for row in refgene_data:
    chrom = row[chr_idx].replace('chr', '')
    if chrom == '23': chrom = 'X'
    if chrom == '24': chrom = 'Y'
    try:
        start = int(row[start_idx])
        end = int(row[end_idx])
        strand = int(row[strand_idx])
    except ValueError:
        continue
    symb = row[symb_idx].strip()
    gene = row[gene_idx].strip()
    if not symb:
        continue
    refgene_normalized.append((chrom, start, end, strand, symb, gene))

# === Load COSMIC ===
with open(cosmic_path, 'r') as f:
    cosmic_lines = [line.strip().split(',') for line in f if line.strip()]
    cosmic_header = cosmic_lines[0]
    gene_col = cosmic_header.index('Gene Symbol')
    cosmic_genes = set(row[gene_col].strip() for row in cosmic_lines[1:] if row[gene_col].strip())

# === Function to annotate top nearby genes ===
def top_nearby_genes(chrom, pos, cosmic_genes, window_cosmic=500_000, window_tight=10_000):
    genes = [row for row in refgene_normalized if row[0] == chrom]
    candidates = []

    for row in genes:
        start, end, symb = row[1], row[2], row[4]
        dist = min(abs(pos - start), abs(pos - end))
        label = 'closest'
        if symb in cosmic_genes and (start >= pos - window_cosmic and end <= pos + window_cosmic):
            label = 'cosmic'
        elif start <= pos + window_tight and end >= pos - window_tight:
            label = 'tight'
        candidates.append((symb, dist, label))

    def sort_key(x):
        priority = {'cosmic': 0, 'tight': 1, 'closest': 2}
        return (priority.get(x[2], 3), x[1])

    candidates_sorted = sorted(candidates, key=sort_key)

    seen = set()
    top_genes = []
    for symb, _, _ in candidates_sorted:
        if symb not in seen:
            top_genes.append(symb)
            seen.add(symb)
        if len(top_genes) == 21:
            break

    while len(top_genes) < 21:
        top_genes.append('')
    return top_genes

# === Prepare output headers ===
gene_i_cols = [f'gene_i_annot_{i}' for i in range(1, 21)]
gene_j_cols = [f'gene_j_annot_{i}' for i in range(1, 21)]
final_header = sig_driv_header + gene_i_cols + ['gene_i_match'] + gene_j_cols + ['gene_j_match']

# === Apply annotation ===
annotated_rows = []
match_i_total = 0
match_j_total = 0
for row in sig_driv_rows:
    row_dict = dict(zip(sig_driv_header, row))

    chrom_i = row_dict['chr_i'].replace('23', 'X').replace('24', 'Y')
    chrom_j = row_dict['chr_j'].replace('23', 'X').replace('24', 'Y')

    try:
        pos_i = int(row_dict['pos_i'])
        pos_j = int(row_dict['pos_j'])
    except ValueError:
        continue

    top5_i = top_nearby_genes(chrom_i, pos_i, cosmic_genes)
    top5_j = top_nearby_genes(chrom_j, pos_j, cosmic_genes)

    gene_i = row_dict.get('gene_i', '').replace('*', '')
    gene_j = row_dict.get('gene_j', '').replace('*', '')
    gene_i_match = 'True' if gene_i in top5_i else 'False'
    gene_j_match = 'True' if gene_j in top5_j else 'False'
 
    if gene_i_match == 'True':
        match_i_total += 1
    if gene_j_match == 'True':
        match_j_total += 1

    final_row = row + top5_i + [gene_i_match] + top5_j + [gene_j_match]
    annotated_rows.append(final_row)

# === Write output ===
with open(output_path, 'w') as out:
    out.write('\t'.join(final_header) + '\n')
    for row in annotated_rows:
        out.write('\t'.join(str(x) for x in row) + '\n')



# === Print % matches ===
percent_match_i = (match_i_total / len(annotated_rows)) * 100
percent_match_j = (match_j_total / len(annotated_rows)) * 100

print(f"Percent exact gene name match (gene_i vs gene_i_annot): {percent_match_i:.2f}%")
print(f"Percent exact gene name match (gene_j vs gene_j_annot): {percent_match_j:.2f}%")
