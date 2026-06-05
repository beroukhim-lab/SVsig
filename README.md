## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/2023.10.13.561748)


## _SVsig_ Install and Dependencies 

SVsig uses MATLAB, which can be obtained [here](https://www.mathworks.com/products/matlab.html). This version has primarily been tested using `MATLAB_R2020a` on macOS Sonoma (14.5) and `MATLAB_R2023b` on x86_64. 

Additionally, install the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

Clone this repo into the directory you wish to run SVSig in. 



## Running _SVsig_


### Preparing files
SVsig accepts input CSV files with required columns identified by **column name** (not fixed position). Additional columns are allowed and ignored unless used downstream for reporting.

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
_SVsig-2Dc_ accounts for novel connections that arise from neiboring rearrangements. 

- To first identify neighboring rearrangements, run [JaBbA](https://github.com/mskilab-org/JaBbA) to obtain a juxtapositions file.

- Open `bin/run2DModel.m`
  - Set the paths to the working directory, rearrangements file you wish to analyze, and output destination file
  - Set the **weights** and **complex** parameters to true. 
  - Run `run2DModel.m`
  - After line 44 in `mix_model_param.m`, run `mix_model_alpha.R`
  - Continue running `mix_model_param.m` until completion

<br>

### Additional parameters (set in run2DModel.m)
- **model_exist**: Boolean to skip model training and use a pre-determined background model. If True, add path to background model in line 23 (complex model) or 25 (simple model) of runSVSig.m. 
- **len_filter**: Only considers rearrangements above this length for calculating significance. Default is 1Mb. 
- **fdr_threshold**: FDR threshold for determining significance. 
- **output_file**: path to output file 
- **complex**: Boolean to run SVSig-2Dc (complex model). 
- **num_breakpoints_per_bin**: Average number of breakpoints within a bin. Determines bin boundaries so that each tile has approximately this number of breakpoints. Currently not used.
- **bin_length**: Length of bin to divide genome. Suggested ranges are 500kb - 2Mb. Note that the number of calculations scales quadratically as bin_length decreases. 
- **weights**: Weight given to each individual connection, ranges from 0-1. Weight=1 for the simple model. For the complex model, weights are obtained from the juxtapositions file after running [JaBbA](https://github.com/mskilab-org/JaBbA)
- **genome_build**: 'hg19' or 'hg_38'.
<br>

### Outputs
_SVsig-2D_ and _SVsig-2Dc_ output a TSV where each row is a rearrangement assigned to a significantly recurrent cluster, denoted with a cluster number. Output columns are:

- **cluster_num**: Recurrent-cluster ID.
- **sid**: Sample ID for the rearrangement.
- **gene_i, gene_j**: Prioritized gene labels near each breakpoint.
  - `**` means the gene pair matches a curated fusion pair reference.
  - `*` means the gene is in the cancer-gene reference list.
  - Genes are ordered by priority: `**` first, then `*`, then unannotated.
- **all_genes_i, all_genes_j**: Full nearby-gene lists before prioritization/truncation.
- **bin_i, bin_j**: Bin IDs for each breakpoint.
- **chr_i, pos_i, strand_i**: Breakpoint i chromosome, position, and strand.
- **chr_j, pos_j, strand_j**: Breakpoint j chromosome, position, and strand.
- **sv_class**: SV class derived from chromosome and strand orientation:
  - `inter_chr`: `chr_i ~= chr_j` (interchromosomal event).
  - `del`: `chr_i == chr_j` and `strand_i = +`, `strand_j = -`.
  - `tandem_dup`: `chr_i == chr_j` and `strand_i = -`, `strand_j = +`.
  - `inv`: `chr_i == chr_j` and `strand_i = strand_j` (`++` or `--`).
  - `unknown`: `chr_i == chr_j` but one or both strands use an unrecognized encoding.
  - Strand parser accepts both symbol and numeric encodings: `+/-`, `1/-1`, and `+1/-1`.
- **tile_qval**: Tile-level BH-FDR q-value.
- **pval**: Tile-level p-value.
- **prob**: Background-model probability for the tile.
- **num_hits**: Number of unique samples in the cluster.
- **stddev_i, stddev_j**: Position standard deviation at each breakpoint across SVs in the cluster.
- **tier**: Tier label based on breakpoint spread thresholds (`tier 1`, `tier 2`, `tier 3`).
<br>

## Tutorial

To ensure that SVsig is installed and running properly, we will run the file `data/TUTORIAL_rearrangements.csv`, which contains a random sampling of 100,000 rearrangements from the dataset in the manuscript. Change the following parameters (use the default for the remaining parameters) and run SVSig-2D:

- **bin_length**: 1e6
- **fdr_threshold**: 0.01

Runtime was measured to be around 7 minutes on a standard laptop with 16GB RAM. The expected output file is shown at `results/TUTORIAL_hitsalljunctions_fdr0.01_1e6bins.txt`. 

To recreate the results in the manuscript from SVSig-2D, use the `data/merged_1.6.1.csv` file, which includes the full set of nearly 300,000 rearrangements from the PCAWG cohort. Additionally, use the following parameters: 

- **bin_length**: 5e5
- **fdr_threshold**: 0.1

## Troubleshooting

Common issues with running _SVSig_ often involve the number of rearrangements in your dataset. SVSig requires a large number of rearrangements since they become sparse once distributed across the genome-wide adjacency matrix. Additionally, at least one rearrangement needs to exist on every chromosome. Ideally, there are at least 100,000 rearrangements in your dataset, although we have run _SVSig_ with data containing only 50,000 rearrangements. For smaller datasets, we recommend increasing the bin_length parameter and increasing the FDR. 

Another option for smaller datsets is generating and loading in the background model using the PCAWG rearrangements (provided in this repo). Afterwards, rearrangements in the new dataset that occur at a higher frequency than the PCAWG background rate can be detected. 



## License
Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute

Please cite: Zhang S, Kumar KH, Shapira O, et al. Detecting significantly recurrent genomic connections from simple and complex rearrangements in the cancer genome. _bioRxiv_ (2023). https://doi.org/10.1101/2023.10.13.561748
