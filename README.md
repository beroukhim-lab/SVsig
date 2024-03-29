## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/2023.10.13.561748)

*MATLAB toolboxes needed:*
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

This version has primarily been tested using `MATLAB_R2020a`. 

## How to run _SVsig_


### Preparing files
SVSig takes in an input .csv file with 15 columns. An example is in this repo under data/merged_1.6.1.csv. Your file must match the column names exactly. 
- **seqnames, start, strand, altchr, altpos, altstrand**: genomic coordinates of both rearrangement breakpoints.
    - Note that chromosome coordinates are integers only. chrX and chrY are changed to 23 and 24, respectively. 
- **dcc_project_code**: histology or tissue type information. 
- **sv_id**: ID for individual rearrangement
- **sid**: Patient ID for the rearrangement
- **donor_unique_id**: Patient ID for the rearrangement
- **weights**: 

Optional columns:
If this information is not available, set column values to arbitrary value. Will not affect ability to run the model.
- **topo**: rearrangement topology information
- **topo_n**: number of rearrangements involved in topology. 
- **mech**: DNA damage repair mechanism predicted to generate the rearrangement. 
- **homseq**: number of base pairs of microhomology at the breakpoint junctions. 
 

### Simple Rearrangements Model (_SVsig-2D_)

SVsig-2D considers each rearrangement to occur independently of each other.
- Open `runSVsig.m`
  - Change path to sample rearrangements file within lines 22-43
- Open `Run2DModel.m` 
  - Make sure: complex, weights, and simulation parameters are all false
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run `Run2DModel.m`

<br>

### Complex Rearrangements Model (_SVsig-2Dc_) 
- Open `runSVsig.m`
  - Change path to sample rearrangements file within lines 22-43
- Open `Run2DModel.m`
  - Set weights and complex parameters to true. 
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run `Run2DModel.m`


There are additional parameters that can be set in Run2DModel.m
- **model_exist**: Boolean to skip model training and use a pre-determined background model. If True, add path to background model in line XX of runSVSig.m
- **len_filter**: Only considers rearrangements above this length for calculating significance. 
- **bks_cluster**:
- **FDR_THRESHOLD**: FDR threshold for determining significance. 
- **output_file**: path to output file 
- **complex**: Boolean to run SVSig-2Dc (complex model). 
- **num_breakpoints_per_bin**: Average number of breakpoints within a bin. To determine bin boundaries. Currently not used. Note to self -- fix this later. 
- **bin_length**: Length of bin to divide genome. Default: 500kb
- **weights**: Weight given to each individual connection, ranges from 0-1. Weight=1 for the simple model. For the complex model, weights are obtained from the juxtapositions file after running [JaBbA](https://github.com/mskilab-org/JaBbA)
- **simulations**: Boolean to test simulated data. 
- **genome_build**: 'hg19' or 'hg_38'.

Note to self: remove local and simulations parameter, add output_file parameter 
ALSO move parameters from break_invasion_model to Run2DModel


### Outputs
_SVsig-2D_ and _SVsig-2Dc_ output a file containing significantly recurrently events. Each unique event is denoted with by a cluster number. The genomic coordinates, subtype, and ID information for each rearrangement in a cluster are displayed. In addition, the following columns are present:  
- **cluster_num**: Cluster number each connection belongs to. 
- **pval**: 
- **prob**: 
- **num_hits**: Number of unique samples containing the rearrangement. 

### Tutorial
make a tutorial later




## License
Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute

Please cite: Zhang S, Kumar KH, Shapira O, et al. Detecting significantly recurrent genomic connections from simple and complex rearrangements in the cancer genome. _bioRxiv_ (2023). https://doi.org/10.1101/2023.10.13.561748
