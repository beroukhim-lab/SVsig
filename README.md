## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/2023.10.13.561748)


## _SVsig_ Install and Dependencies 

SVsig uses MATLAB, which can be obtained [here](https://www.mathworks.com/products/matlab.html). This version has primarily been tested using `MATLAB_R2020a` on macOS Sonoma (14.5) and `MATLAB_R2023b` on x86_64. 

Additionally, install the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

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

No additional columns are required. Inputs may contain extra columns, and SVsig will ignore them for model fitting.

<br>

### Simple Rearrangements Model (_SVsig-2D_)

_SVsig-2D_ considers each rearrangement to occur independently of each other.
- Open `bin/run2DModel.m`
  - Set the paths to the working directory, rearrangements file you wish to analyze, and output destination file
  - Ensure: **complex**, **weights**, and **model_exist** parameters are false
  - Run `run2DModel.m`

<br>

### Complex Rearrangements Model (_SVsig-2Dc_) 
_SVsig-2Dc_ accounts for novel connections that arise from neighboring rearrangements. 

- To first identify neighboring rearrangements, run [JaBbA](https://github.com/mskilab-org/JaBbA) to obtain a juxtapositions file.

- Open `bin/run2DModel.m`
  - Set the paths to the working directory, rearrangements file you wish to analyze, and output destination file
  - Set the **weights** and **complex** parameters to true. 
  - Run `run2DModel.m`
  - After line 44 in `mix_model_param.m`, run `mix_model_alpha.R`
  - Continue running `mix_model_param.m` until completion

<br>

### Additional parameters (set in run2DModel.m)
- **model_exist**: Boolean to skip model training and use a pre-determined background model. If True, set the **background_model_path** parameter to point to a pre-trained model file (for either simple or complex model mode). 
- **len_filter**: Only considers rearrangements above this length for calculating significance. Default is 1Mb. 
- **fdr_threshold**: FDR threshold for determining significance. 
- **output_file**: Path to output file. 
- **complex**: Boolean to run _SVsig-2Dc_ (complex model). 
- **num_breakpoints_per_bin**: Average number of breakpoints within a bin. Determines bin boundaries so that each tile has approximately this number of breakpoints. Currently not used.
- **bin_length**: Length of bin to divide genome. Suggested ranges are 500kb - 2Mb. Note that the number of calculations scales quadratically as bin_length decreases. 
- **weights**: Weight given to each individual connection, ranges from 0-1. Weight=1 for the simple model. For the complex model, weights are obtained from the juxtapositions file after running [JaBbA](https://github.com/mskilab-org/JaBbA).
- **genome_build**: 'hg19' or 'hg38'.
- **n_binshifts**: Number of additional bin-shift runs to perform. Default is 0 (no bin-shifting). Each shift tests the robustness of cluster detection by offsetting bin boundaries. When n_binshifts > 0, results from all runs are merged into superclusters based on spatial proximity.
- **max_workers_override**: Maximum number of parallel workers for bin-shift runs. If empty or 0, automatically detected from system CPU count. Only used when n_binshifts > 0.
- **supercluster_distance_threshold**: Distance threshold (in bp) for merging clusters across bin-shifted runs. Default is 50000 bp. Clusters are considered part of the same supercluster if both breakpoints are within this distance on both breakpoint coordinates.
<br>

### Bin-Shifting Mode

When `n_binshifts > 0`, _SVsig_ runs the model multiple times with different bin boundary offsets to assess cluster robustness:
- **Run 0** (baseline): Standard binning with no offset.
- **Runs 1 to n_binshifts**: Each run applies a fractional offset k/(n_binshifts+1) of the bin_length, creating alternative binning schemes.
- **Merging**: After filtering per-run results, clusters from all runs are merged into superclusters based on proximity thresholds on both breakpoint coordinates.
- **Output**: Single merged table with supercluster IDs, per-SV cluster tracking (sv_cluster_ids, sv_cluster_tiers), and multi-run statistics.
_SVsig-2D_ and _SVsig-2Dc_ output a TSV where each row is a rearrangement assigned to a significantly recurrent cluster. Output columns are:

**Cluster and Supercluster Identification:**
- **supercluster_id**: Integer ID for the merged cluster (same run or across bin-shifted runs). In single-run mode, this is a remapped cluster ID; in multi-run mode, it represents the connected component ID after merging.
- **sv_cluster_ids**: Tuple format identifying all clusters this SV belongs to, e.g., `(0, 5), (1, 3)`. Format is `(run_idx, cluster_num)` for each cluster.
- **sv_cluster_tiers**: Tier values for each cluster listed in `sv_cluster_ids`, in matching order (e.g., `tier 1, tier 2`). Useful for tracking cluster precision before merging.
- **supercluster_cluster_ids**: All clusters in this supercluster in tuple format, e.g., `(0, 3), (0, 5), (1, 2)`. Shows all (run_idx, cluster_num) pairs in the supercluster.
- **supercluster_tier**: Tier label for the supercluster based on position spread thresholds: `tier 1` (both stddev ≤ cutoff), `tier 2` (one stddev ≤ cutoff), `tier 3` (both stddev > cutoff).
- **shifts_found_count**: Number of bin-shift runs where this SV was detected. Value is 1 for single-run (no bin-shifting), ≥1 for multi-run.

**Sample and Gene Information:**
- **sid**: Sample ID for the rearrangement.
- **gene_i, gene_j**: Prioritized gene labels near each breakpoint.
  - `**` means the gene pair matches a curated fusion pair reference.
  - `*` means the gene is in the cancer-gene reference list.
  - Genes are ordered by priority: `**` first, then `*`, then unannotated.
- **all_genes_i, all_genes_j**: Full nearby-gene lists before prioritization/truncation.

**Genomic Coordinates:**
- **bin_i, bin_j**: Bin IDs for each breakpoint.
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
- **tile_qval, pval, prob**: Per-run tile q-value, p-value, and background probability (from the specific run/cluster where SV was detected).
- **sv_tile_qval, sv_pval, sv_prob**: Comma-separated lists (multi-run only) matching `sv_cluster_ids` order; per-cluster statistics for each run where the SV appears.
- **supercluster_pval**: Minimum p-value across all clusters in supercluster.
- **supercluster_prob**: Mean probability across all SVs in supercluster.
- **supercluster_nsamples**: Number of unique samples in supercluster.
- **supercluster_tier**: Tier label (`tier 1/2/3`) based on position spread.
- **stddev_i, stddev_j**: Position standard deviation at each breakpoint across supercluster SVs.
<br>

## Tutorial

To ensure that _SVsig_ is installed and running properly, we will run the file `data/TUTORIAL_rearrangements.csv`, which contains a random sampling of 100,000 rearrangements from the dataset in the manuscript. Change the following parameters (use the default for the remaining parameters) and run _SVsig_-2D:

- **bin_length**: 1e6
- **fdr_threshold**: 0.01

Runtime was measured to be around 7 minutes on a standard laptop with 16GB RAM. The expected output file is shown at `results/TUTORIAL_hitsalljunctions_fdr0.01_1e6bins.txt`. 

To recreate the results in the manuscript from _SVsig-2D_, use the `data/merged_1.6.1.csv` file, which includes the full set of nearly 300,000 rearrangements from the PCAWG cohort. Additionally, use the following parameters: 

- **bin_length**: 5e5
- **fdr_threshold**: 0.1

## Troubleshooting

Common issues with running _SVsig_ often involve the number of rearrangements in your dataset. _SVsig_ requires a large number of rearrangements since they become sparse once distributed across the genome-wide adjacency matrix. Additionally, at least one rearrangement needs to exist on every chromosome. Ideally, there are at least 100,000 rearrangements in your dataset, although we have run _SVsig_ with data containing only 50,000 rearrangements. For smaller datasets, we recommend increasing the bin_length parameter and increasing the FDR. 

Another option for smaller datasets is generating and loading in the background model using the PCAWG rearrangements (provided in this repo). Afterwards, rearrangements in the new dataset that occur at a higher frequency than the PCAWG background rate can be detected. 



## License
Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute

Please cite: Zhang S, Kumar KH, Shapira O, et al. Detecting significantly recurrent genomic connections from simple and complex rearrangements in the cancer genome. _bioRxiv_ (2023). https://doi.org/10.1101/2023.10.13.561748
