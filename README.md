## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/2023.10.13.561748)


## _SVsig_ Install and Dependencies 

SVsig uses MATLAB, which can be obtained [here](https://www.mathworks.com/products/matlab.html). This version has primarily been tested using `MATLAB_R2020a` on macOS Sonoma (14.5). 

Additionally, install the following MATLAB toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

Clone this repo into the directory you wish to run SVSig in. 



## Running _SVsig_


### Preparing files
SVsig takes in an input .csv file with 15 columns. An example is in this repo under `data/merged_1.6.1.csv`. Your file must match the column names exactly. 
- **seqnames, start, strand, altchr, altpos, altstrand**: genomic coordinates of both rearrangement breakpoints.
    - Note that chromosome coordinates are integers only. chrX and chrY are changed to 23 and 24, respectively. 
- **dcc_project_code**: histology or tissue type information. 
- **sv_id**: ID for individual rearrangement
- **sid**: Patient ID for the rearrangement
- **donor_unique_id**: Patient ID for the rearrangement
- **weights**: Ranges from 0-1, representing the weight an individual connection is given to an entire rearragnement event. 

Optional columns:
If this information is not available, set column values to arbitrary value. Will not affect ability to run the model.
- **topo**: rearrangement topology information
- **topo_n**: number of rearrangements involved in topology. 
- **mech**: DNA damage repair mechanism predicted to generate the rearrangement. 
- **homseq**: number of base pairs of microhomology at the breakpoint junctions. 

<br>

### Simple Rearrangements Model (_SVsig-2D_)

_SVsig-2D_ considers each rearrangement to occur independently of each other.
- Open `runSVsig.m`
  - Change path to sample rearrangements file within lines 22-43
- Open `Run2DModel.m` 
  - Make sure: complex, weights, and simulation parameters are all false
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run `Run2DModel.m`

<br>

### Complex Rearrangements Model (_SVsig-2Dc_) 
_SVsig-2Dc_ accounts for novel connections that arise from neiboring rearrangements. 

- To first -----, run [JaBbA](https://github.com/mskilab-org/JaBbA), which 


- Open `runSVsig.m`
  - Change path to sample rearrangements file within lines 22-43
- Open `Run2DModel.m`
  - Set weights and complex parameters to true. 
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run `Run2DModel.m`

<br>

### Additional parameters (set in Run2DModel.m)
- **model_exist**: Boolean to skip model training and use a pre-determined background model. If True, add path to background model in line 23 (complex model) or 25 (simple model) of runSVSig.m. 
- **len_filter**: Only considers rearrangements above this length for calculating significance. Default is 1Mb. 
- **bks_cluster**: Set to 1. 
- **FDR_THRESHOLD**: FDR threshold for determining significance. 
- **output_file**: path to output file 
- **complex**: Boolean to run SVSig-2Dc (complex model). 
- **num_breakpoints_per_bin**: Average number of breakpoints within a bin. Determines bin boundaries so that each tile has approximately this number of breakpoints. Currently not used.
- **bin_length**: Length of bin to divide genome. Suggested ranges are 500kb - 2Mb. Note that the number of calculations scales quadratically as bin_length decreases. 
- **weights**: Weight given to each individual connection, ranges from 0-1. Weight=1 for the simple model. For the complex model, weights are obtained from the juxtapositions file after running [JaBbA](https://github.com/mskilab-org/JaBbA)
- **simulations**: Boolean to test simulated data. 
- **genome_build**: 'hg19' or 'hg_38'.


### Outputs
_SVsig-2D_ and _SVsig-2Dc_ output a file containing significantly recurrently events. Each unique event is denoted with by a cluster number. The genomic coordinates, subtype, and ID information for each rearrangement in a cluster are displayed. In addition, the following columns are present:  
- **cluster_num**: Cluster number each connection belongs to. 
- **pval**: Significance for the rearrangement event. 
- **prob**: 
- **num_hits**: Number of unique samples containing the rearrangement. 

## Tutorial

To ensure that SVsig is installed and running properly, we will run the file in `data/TUTORIAL_rearrangements.csv`. Change the following parameters:

- **bin_length**: 5e6
- **FDR_THRESHOLD**: 0.01

Runtime was measured to be around 7 minutes on a standard laptop with 16GB RAM. The expected output file is shown at `results/TUTORIAL_hitsalljunctions_fdr0.01_1e6bins.txt`. 

To recreate the results in the manuscript from SVSig-2D, use the `data/merged_1.6.1.csv` file, which includes the 300,000 rearrangements from the PCAWG cohort. Additionally, use the following parameters: 

- **bin_length**: 5e5
- **FDR_THRESHOLD**: 0.1

## Troubleshooting

Common issues with running _SVSig_ often involve the number of rearrangements in your dataset. SVSig requires a large number of rearrangements since they become sparse once distributed across the genome-wide adjacency matrix. Additionally, at least one rearrangement needs to exist on every chromosome. Ideally, there are at least 100,000 rearrangements in your dataset, although we have run _SVSig_ with data containing only 50,000 rearrangements. For smaller datasets, we recommend increasing the bin_length parameter and increasing the FDR. 

Another option for smaller datsets is generating and loading in the background model using the PCAWG rearrangements (provided in this repo). Afterwards, rearrangements in the new dataset that occur at a higher frequency than the PCAWG background rate can be detected. 



## License
Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute

Please cite: Zhang S, Kumar KH, Shapira O, et al. Detecting significantly recurrent genomic connections from simple and complex rearrangements in the cancer genome. _bioRxiv_ (2023). https://doi.org/10.1101/2023.10.13.561748
