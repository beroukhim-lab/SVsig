## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/2023.10.13.561748)


## _SVsig_ Install and Dependencies 

_SVsig_ uses MATLAB, which can be obtained [here](https://www.mathworks.com/products/matlab.html). This version has primarily been tested using `MATLAB_R2020a` on macOS Sonoma (14.5), `MATLAB_R2023b` on x86_64, and `MATLAB_R2024b` on macOS.

Additionally, install the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox (optional; used to run bin-shifted runs in parallel, otherwise they run serially)

Clone this repo into the directory you wish to run _SVsig_ in. 



## Running _SVsig_


### Preparing files
_SVsig_ accepts input CSV files with required columns identified by **column name**. Additional columns are allowed but ignored.

Required columns (both simple and complex models):
- **seqnames, start, strand, altchr, altpos, altstrand**: genomic coordinates and strands of both rearrangement breakpoints.
  - `seqnames` and `altchr` can be `1-22`, `X`, `Y`, `chr1-22`, `chrX`, or `chrY`.
- **sid**: sample identifier used for sample-level de-duplication and per-cluster recurrence counts.

Conditionally required:
- **weights**: required for meaningful complex-model runs (`complex_model=true`).
  - For simple model runs, if `weights` is missing, SVsig defaults all weights to `1`.

No additional columns are required. Inputs may contain extra columns, and _SVsig_ will ignore them for model fitting.

<br>

### Run With cfg (MATLAB)

`run2DModel` accepts a configuration struct, so you do not need to edit source files.

Simple model (_SVsig-2D_):

```matlab
cfg = struct();
cfg.work_dir = pwd;
cfg.sv_file = 'path/to/input.csv';
cfg.output_file = 'path/to/output.tsv';

cfg.complex_model = false;
cfg.weights = false;
cfg.model_exist = false;

hits = run2DModel(cfg);
```

Complex model (_SVsig-2Dc_):

```matlab
cfg = struct();
cfg.work_dir = pwd;
cfg.sv_file = 'path/to/input.csv';
cfg.output_file = 'path/to/output.tsv';

cfg.complex_model = true;
cfg.weights = true;
cfg.model_exist = false;

hits = run2DModel(cfg);
```

If you already have a precomputed model:

```matlab
cfg.model_exist = true;
cfg.model_file = 'path/to/background_model.mat';
```

The output delimiter is chosen from the `output_file` extension: `.csv` writes a comma-delimited file, and `.tsv` (or any other extension) writes a tab-delimited file. For example:

```matlab
cfg.output_file = 'path/to/output.csv';  % comma-delimited
```

<br>

### Run With CLI Flags

`run2DModel` also supports flags (useful for compiled/batch runs). For example:

```bash
run2DModel -sv path/to/input.csv -out path/to/output.tsv -bin 1e6 -fdr 0.01 -bix true
```

Common flags:
- `-sv`, `--sv_file` input CSV (required)
- `-out`, `--output_file` output TSV (required)
- `-wd`, `--work_dir` project root
- `-complex`, `--complex_model` true/false
- `-weights`, `--weights` true/false
- `-model_exist`, `--model_exist` true/false
- `-model`, `--model_file` precomputed model path
- `-bin`, `--bin_length` bin length
- `-fdr`, `--fdr_threshold` BH-FDR cutoff
- `-nshift`, `--n_binshifts` number of additional shifted runs
- `-maxw`, `--max_workers` parallel worker cap (`1` forces serial execution)
- `-scdist`, `--supercluster_distance_threshold` supercluster merge distance (bp)
- `-genome`, `--genome_build` genome build (`hg_19` or `hg_38`; selects `data/tracks/<genome_build>/`)
- `-h`, `--help` print full option help

<br>

### Hyperparameters (cfg keys / CLI equivalents)

Core model:
- `complex_model` (`-complex`): enable complex model mode
- `weights` (`-weights`): use per-rearrangement weights
- `model_exist` (`-model_exist`): load precomputed background model
- `model_file` (`-model`): background model path when `model_exist=true`

Filtering and significance:
- `fdr_threshold` (`-fdr`): BH-FDR significance cutoff
- `len_filter` (`-len`): minimum SV span used in filtering
- `std_filter` (`-std`): recurrence spread filter threshold
- `tier_std_cutoff` (`-tier`): cutoff used for tier assignment
- `low_density_threshold` (`-ldt`): low-density bin removal threshold

Binning and superclustering:
- `bin_length` (`-bin`): bin size in bp (default `5e5`)
- `num_breakpoints_per_bin` (`-nbpb`): retained for compatibility but currently unused; bins are fixed-width and set by `bin_length`
- `n_binshifts` (`-nshift`): additional shifted runs (default `0`)
- `supercluster_distance_threshold` (`-scdist`): max breakpoint distance for merging clusters

Runtime and tracks:
- `max_workers` (`-maxw`): max parallel workers
- `chr_list` (`-chr`): chromosomes to include
- `genome_build` (`-genome`): genome build for track resolution
- `random_seed` (`-seed`): RNG seed. If unset, every run starts from MATLAB's default seed (`rng('default')`). If set, run `k` (0-based) uses `seed + k`
- `save_bin_index` (`-bix`): write bin-index sidecar files. Names insert a tag before the output extension, e.g. `out.tsv` gives `out.bin_indices.tsv` (single run) or `out.bin_indices.shift<N>_bp<S>.tsv` (per shift run)
- `save_per_run_output` (`-per_run`): when `n_binshifts > 0`, also write one output file per shift run (e.g. `out.tsv` gives `out.shift<N>_bp<S>.tsv`)
- `ref_genes_file`, `cancer_genes_file`, `curated_fusions_file`, `chrom_sizes_file`, `blacklist_file`, `l1_elements_file`: optional track file overrides

<br>

### Bin-Shifting Mode

When `n_binshifts > 0`, _SVsig_ runs the model multiple times with different bin boundary offsets to assess cluster robustness:
- **Run 0** (baseline): Standard binning with no offset.
- **Runs 1 to n_binshifts**: Each run translates the bin lattice by `floor(k/(n_binshifts+1) * bin_length)` bp. Because the translated lattice would otherwise leave positions `1..shift` uncovered, each chromosome also gets a leading partial bin spanning `1` to `shift`.
- **Execution**: Runs execute in parallel with `parfor` when the Parallel Computing Toolbox is available and more than one worker is allowed (workers = min(runs, CPU cores, `max_workers`)); otherwise, or if the parallel run fails, they run serially.
- **Merging**: After filtering per-run results, clusters from all runs are merged into superclusters based on proximity thresholds on both breakpoint coordinates.
- **Output**: Single merged table with one row per unique SV, supercluster IDs, per-SV cluster tracking (`sv_cluster_ids`, `sv_cluster_tiers`), and multi-run statistics. Per-run artifacts (`save_per_run_output`, `save_bin_index`) are written as each run completes.

### Outputs

_SVsig-2D_ and _SVsig-2Dc_ output a table (comma-delimited if `output_file` ends in `.csv`, otherwise tab-delimited) where each row is a rearrangement assigned to a significantly recurrent cluster. Output columns are:

**Cluster and Supercluster Identification:**
- **supercluster_id**: Integer ID for the merged cluster (same run or across bin-shifted runs). In single-run mode, this is a remapped cluster ID; in multi-run mode, it represents the connected component ID after merging.
- **sv_cluster_ids**: Tuple format identifying all clusters this SV belongs to, e.g., `(0, 5), (1, 3)`. Format is `(run_idx, cluster_num)` for each cluster.
- **sv_cluster_tiers**: Tier values for each cluster listed in `sv_cluster_ids`, in matching order (e.g., `tier 1, tier 2`). Useful for tracking cluster precision before merging.
- **supercluster_cluster_ids**: All clusters in this supercluster in tuple format, e.g., `(0, 3), (0, 5), (1, 2)`. Shows all (run_idx, cluster_num) pairs in the supercluster.
- **supercluster_tier**: Tier label for the supercluster based on position spread thresholds: `tier 1` (both stddev ≤ cutoff), `tier 2` (one stddev ≤ cutoff), `tier 3` (both stddev > cutoff).
- **shifts_found_count**: Number of bin-shift runs where this SV was detected. This is a supercluster-level count (the same for every SV in a supercluster), not a per-SV count. Value is 1 for single-run (no bin-shifting).

**Sample and Gene Information:**
- **sid**: Sample ID for the rearrangement.
- **gene_i, gene_j**: Prioritized gene labels near each breakpoint.
  - `**` means the gene pair matches a curated fusion pair reference.
  - `*` means the gene is in the cancer-gene reference list.
  - Genes are ordered by priority: `**` first, then `*`, then unannotated.
- **all_genes_i, all_genes_j**: Full nearby-gene lists before prioritization/truncation.

**Genomic Coordinates:**
- **sv_bin_i, sv_bin_j**: Comma-separated bin IDs for each cluster where the SV appears, in matching order with `sv_cluster_ids`.
- **chr_i, pos_i, strand_i**: Breakpoint i chromosome, position, and strand.
- **chr_j, pos_j, strand_j**: Breakpoint j chromosome, position, and strand.
- **sv_class**: SV class derived from chromosome and strand orientation:
  - `inter_chr`: `chr_i ≠ chr_j` (interchromosomal event).
  - `del`: `chr_i == chr_j` and `strand_i = +`, `strand_j = -`.
  - `tandem_dup`: `chr_i == chr_j` and `strand_i = -`, `strand_j = +`.
  - `inv`: `chr_i == chr_j` and `strand_i = strand_j` (`++` or `--`).
  - `unknown`: `chr_i == chr_j` but one or both strands use an unrecognized encoding.
  - Strand parser accepts both symbol and numeric encodings: `+/-`, `1/-1`, and `+1/-1`.

**Statistical Measures:**
- **tile_qval, pval, prob**: In single-run mode, the tile q-value, p-value, and background probability of the cluster. In multi-run mode these are supercluster-level: minimum q-value, minimum p-value, and mean probability across the supercluster.
- **sv_tile_qval, sv_pval, sv_prob**: Comma-separated lists matching `sv_cluster_ids` order; per-cluster statistics for each run where the SV appears (a single value in single-run mode).
- **supercluster_pval**: Minimum p-value across all clusters in supercluster.
- **supercluster_prob**: Mean probability across all SVs in supercluster.
- **supercluster_nsamples**: Number of unique samples in supercluster.
- **supercluster_tier**: Tier label (`tier 1/2/3`) based on position spread.
- **stddev_i, stddev_j**: Position standard deviation at each breakpoint across supercluster SVs.

<br>

## Reference Genome Files and Formats

_SVsig_ requires six reference files for gene annotation, fusion detection, and masking. These files are resolved automatically from `data/tracks/<genome_build>/` based on your `genome_build` parameter, or can be explicitly overridden via cfg or CLI flags.

### File Format Specifications

**Chromosome sizes:**
- Plain text, tab- or space-separated.
- Two columns: chromosome identifier and length (bp).
- Headerless files are accepted; optional header row will be skipped if present.
- Chromosome names: accepts numeric (1-24), string (chr1-22, chrX, chrY), or mixed formats.
- Example (headerless): `1  249250621` or `chrX  155270560`

**Reference genome:**
- Tab-delimited text with required header row (.txt or .txt.gz).
- Required columns (flexible names):
  - Chromosome: `chr`, `chrom`, or `chromosome`
  - Start: `start`, `txStart`, or `chromStart`
  - End: `end`, `txEnd`, `chromEnd`, or `xEnd`
  - Gene symbol: `gene_symbol`, `name2`, `gene`, or `symbol`
- Additional columns are ignored.
- Example header: `chr	start	end	gene_symbol`

**Cancer genes:**
- Tab- or comma-delimited text (.tsv, .csv, or .txt); can be gzipped.
- Single gene-name column, with optional header.
- Flexible column names: looks for `gene_symbol`, `gene`, `symbol`, or `Gene Symbol`.
- If no matching column header is found, uses the first column.
- Fallback: plain text file (one gene name per line; lines starting with `#` are ignored).
- Example with header: `gene_symbol` followed by gene names.
- Example headerless: `TP53` on line 1, `BRCA1` on line 2, etc.

**Curated fusions:**
- Tab-delimited text with required header row (.tsv or .txt, can be gzipped).
- Two columns for gene pairs with flexible names:
  - First gene: `gene_A`, `gene1`, or `5prime_gene`
  - Second gene: `gene_B`, `gene2`, or `3prime_gene`
- Example header: `gene_A	gene_B`

**Blacklist regions:**
- BED format (.bed or .bed.gz): tab-separated columns [chromosome, start, end, ...].
- Headerless (BED standard).
- Coordinates are 0-based, half-open [start, end) per BED specification.
- Example: `1	1000	2000` or `chr1	1000	2000`

**L1 elements:**
- Same format as blacklist regions (BED3 or BED+, headerless, 0-based coordinates).
- Example: `chrX	500	1500`

### Chromosome Naming

All file loaders accept mixed chromosome naming conventions and automatically normalize them:
- Numeric: `1`, `2`, `23` (X), `24` (Y)
- String with chr prefix: `chr1`, `chr2`, `chrX`, `chrY`
- String without prefix: `1`, `2`, `X`, `Y`
- Case-insensitive: `X`, `x`, `CHR23` all map to chromosome 23.

Only canonical chromosomes 1-22, X (23), and Y (24) are kept; other entries are ignored.

### Default Sources

The repository currently ships default reference files for hg19 only (`data/tracks/hg_19/`). For hg38, set `genome_build = 'hg_38'` and either place files in `data/tracks/hg_38/` using the default filenames in `bin/resolve_track_paths.m`, or pass them via the override options below.

*TODO: Haruna needs to fill this out with sources and citations!* 

| Input File | hg19 (GRCh37) | hg38 (GRCh38) |
|---|---|---|
| Chromosome sizes (`chrom_sizes_file`) | **File:** `data/tracks/hg_19/hg19.chrom.sizes`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |
| Reference genome (`ref_genes_file`) | **File:** `data/tracks/hg_19/hg19.ncbiRefSeq.colnames.txt`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |
| Cancer genes (`cancer_genes_file`) | **File:** `data/tracks/hg_19/Cosmic_CancerGeneCensus_v104_GRCh37.tsv`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |
| Curated fusions (`curated_fusions_file`) | **File:** `data/tracks/hg_19/CuratedFusionGene.colnames.txt`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |
| Blacklist regions (`blacklist_file`) | **File:** `data/tracks/hg_19/mask_track.bed`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |
| L1 elements (`l1_elements_file`) | **File:** `data/tracks/hg_19/l1_track.bed`<br>**Source/Citation:** _[ADD LINK/CITATION]_ | **File:** `data/tracks/hg_38/[ADD_FILENAME]`<br>**Source/Citation:** _[ADD LINK/CITATION]_ |

### Overriding Track Files

To use custom or alternative track files, specify them via cfg or CLI:

**MATLAB cfg:**
```matlab
cfg.ref_genes_file = '/path/to/custom_refgenes.txt';
cfg.cancer_genes_file = '/path/to/custom_cancer_genes.tsv';
cfg.curated_fusions_file = '/path/to/custom_fusions.tsv';
cfg.chrom_sizes_file = '/path/to/custom.chrom.sizes';
cfg.blacklist_file = '/path/to/custom_blacklist.bed';
cfg.l1_elements_file = '/path/to/custom_l1.bed';
```

**CLI flags:**
```bash
run2DModel -sv input.csv -out output.tsv \
  -ref_genes /path/to/custom_refgenes.txt \
  -cancer_genes /path/to/custom_cancer_genes.tsv \
  -curated_fusions /path/to/custom_fusions.tsv \
  -chrom_sizes /path/to/custom.chrom.sizes \
  -blacklist /path/to/custom_blacklist.bed \
  -l1 /path/to/custom_l1.bed
```

If any override is not specified, _SVsig_ falls back to the default path for the specified `genome_build`.

<br>

## Tutorial

To verify your installation, run _SVsig-2D_ on `data/TUTORIAL_rearrangements.csv` (100,000 sampled rearrangements):

```matlab
cfg = struct();
cfg.work_dir = pwd;
cfg.sv_file = 'data/TUTORIAL_rearrangements.csv';
cfg.output_file = 'results/tutorial_output.tsv';
cfg.complex_model = false;
cfg.model_exist = false;
cfg.bin_length = 1e6;
cfg.fdr_threshold = 0.01;
run2DModel(cfg);
```

Equivalent CLI invocation:

```bash
run2DModel -sv data/TUTORIAL_rearrangements.csv -out results/tutorial_output.tsv -complex false -model_exist false -bin 1e6 -fdr 0.01
```

Runtime is typically around 7 minutes on a standard laptop with 16GB RAM. A representative expected output is available at `results/TUTORIAL_hitsalljunctions_fdr0.01_1e6bins.txt`.

To recreate manuscript-scale _SVsig-2D_ runs with `data/merged_1.6.1.csv` (~300,000 rearrangements), use:

```matlab
cfg = struct();
cfg.work_dir = pwd;
cfg.sv_file = 'data/merged_1.6.1.csv';
cfg.output_file = 'results/merged_1.6.1_hits.tsv';
cfg.complex_model = false;
cfg.model_exist = false;
cfg.bin_length = 5e5;
cfg.fdr_threshold = 0.1;
run2DModel(cfg);
```

## Troubleshooting

Common issues with running _SVsig_ often involve the number of rearrangements in your dataset. _SVsig_ performs best with large cohorts of rearrangements because events become sparse once distributed across the genome-wide adjacency matrix. Ideally, there are at least 100,000 rearrangements in your dataset, although we have run _SVsig_ with data containing only 50,000 rearrangements. For smaller datasets, we recommend increasing the bin_length parameter and increasing the FDR. 

Another option for smaller datasets is generating and loading in the background model using the PCAWG rearrangements (provided in this repo) or another pre-generated background model complementary to your dataset. Afterwards, rearrangements in the new dataset that occur at a higher frequency than the PCAWG (or other) background rate can be detected. 



## License
Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute

Please cite: Zhang S, Kumar KH, Shapira O, et al. Detecting significantly recurrent genomic connections from simple and complex rearrangements in the cancer genome. _bioRxiv_ (2023). https://doi.org/10.1101/2023.10.13.561748
