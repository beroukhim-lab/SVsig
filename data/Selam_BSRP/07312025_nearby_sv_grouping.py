import pandas as pd
from pyliftover import LiftOver
import argparse

def lift_coord(lo, chrom, pos):
    try:
        pos_int = int(float(pos))  # robust int conversion
    except:
        return None, None
    result = lo.convert_coordinate('chr' + str(chrom), pos_int)
    if result:
        new_chr, new_pos, _, _ = result[0]
        clean_chr = new_chr.replace('chr', '')
        return clean_chr, int(new_pos)
    return None, None

def liftover_pcawg(pcawg_df, chain_file):
    lo = LiftOver(chain_file)
    lifted_i = pcawg_df.apply(lambda row: lift_coord(lo, row['chr_i'], row['pos_i']), axis=1, result_type='expand')
    pcawg_df['chr_i_hg38'] = lifted_i[0]
    pcawg_df['pos_i_hg38'] = lifted_i[1]

    lifted_j = pcawg_df.apply(lambda row: lift_coord(lo, row['chr_j'], row['pos_j']), axis=1, result_type='expand')
    pcawg_df['chr_j_hg38'] = lifted_j[0]
    pcawg_df['pos_j_hg38'] = lifted_j[1]

    # drop rows with liftover failure
    pcawg_df = pcawg_df.dropna(subset=['pos_i_hg38', 'pos_j_hg38'])
    # ensure types
    pcawg_df['pos_i_hg38'] = pcawg_df['pos_i_hg38'].astype(int)
    pcawg_df['pos_j_hg38'] = pcawg_df['pos_j_hg38'].astype(int)
    pcawg_df['chr_i_hg38'] = pcawg_df['chr_i_hg38'].astype(str)
    pcawg_df['chr_j_hg38'] = pcawg_df['chr_j_hg38'].astype(str)

    return pcawg_df

def main():
    parser = argparse.ArgumentParser(description="Compare CCLE and PCAWG SV breakpoints after liftover.")
    parser.add_argument('--ccle', required=True, help='CCLE breakpoints file (hg38)')
    parser.add_argument('--pcawg', required=True, help='PCAWG breakpoints file (hg19)')
    parser.add_argument('--chain', required=True, help='hg19 to hg38 chain file')
    parser.add_argument('--sep', default='\t', help='Delimiter for input files (default tab)')
    args = parser.parse_args()

    # Load CCLE (already hg38)
    driv = pd.read_csv(args.ccle, sep=args.sep)
    # Load PCAWG (hg19, needs liftover)
    pcawg = pd.read_csv(args.pcawg, sep=args.sep)

    # Liftover PCAWG from hg19 to hg38
    pcawg = liftover_pcawg(pcawg, args.chain)

    # Ensure CCLE chr columns are strings
    driv['chr_i'] = driv['chr_i'].astype(str)
    driv['chr_j'] = driv['chr_j'].astype(str)

    same_bin_results = []

    for _, d_row in driv.iterrows():
        chr_i_ccle = d_row['chr_i']
        chr_j_ccle = d_row['chr_j']
        pos_i_ccle = d_row['pos_i']
        pos_j_ccle = d_row['pos_j']
        gene_i_ccle = d_row['gene_i']
        gene_j_ccle = d_row['gene_j']

        for _, s_row in pcawg.iterrows():
            chr_i_pcawg = s_row['chr_i_hg38']
            chr_j_pcawg = s_row['chr_j_hg38']
            pos_i_pcawg = s_row['pos_i_hg38']
            pos_j_pcawg = s_row['pos_j_hg38']
            gene_i_pcawg = s_row['gene_i']
            gene_j_pcawg = s_row['gene_j']

            # Match chromosomes either direct or swapped
            if chr_i_ccle == chr_i_pcawg and chr_j_ccle == chr_j_pcawg:
                d1 = abs(pos_i_ccle - pos_i_pcawg)
                d2 = abs(pos_i_ccle - pos_j_pcawg)
                if d1 <= d2:
                    min_i_dist = d1
                    dist_j = abs(pos_j_ccle - pos_j_pcawg)
                else:
                    min_i_dist = d2
                    dist_j = abs(pos_j_ccle - pos_i_pcawg)
            elif chr_i_ccle == chr_j_pcawg and chr_j_ccle == chr_i_pcawg:
                d1 = abs(pos_i_ccle - pos_j_pcawg)
                d2 = abs(pos_i_ccle - pos_i_pcawg)
                if d1 <= d2:
                    min_i_dist = d1
                    dist_j = abs(pos_j_ccle - pos_i_pcawg)
                else:
                    min_i_dist = d2
                    dist_j = abs(pos_j_ccle - pos_j_pcawg)
            else:
                continue

            total_dist = min_i_dist + dist_j
            in_same_bin = total_dist <= 2_000_000  # 2 Mb threshold

            same_bin_results.append({
                'chr_ccle_i': chr_i_ccle,
                'chr_ccle_j': chr_j_ccle,
                'ccle_gene_i': gene_i_ccle,
                'ccle_gene_j': gene_j_ccle,
                'ccle_pos_i': pos_i_ccle,
                'ccle_pos_j': pos_j_ccle,
                'chr_pcawg_i': chr_i_pcawg,
                'chr_pcawg_j': chr_j_pcawg,
                'pcawg_gene_i': gene_i_pcawg,
                'pcawg_gene_j': gene_j_pcawg,
                'pcawg_pos_i': pos_i_pcawg,
                'pcawg_pos_j': pos_j_pcawg,
                'subtype_pcawg': s_row.get('subtype', ''),
                'subtype_ccle': d_row.get('subtype', ''),
                'min_i_dist': min_i_dist,
                'dist_j': dist_j,
                'total_dist': total_dist,
                'same_bin': in_same_bin
            })

    bin_df = pd.DataFrame(same_bin_results)

    print(f"bins: {len(bin_df[bin_df['same_bin']])}")
    print(f"not_bins: {len(bin_df[~bin_df['same_bin']])}")

    # Optionally save or display matched bins
    matched = bin_df[bin_df['same_bin']]
    print(matched)

if __name__ == '__main__':
    main()

