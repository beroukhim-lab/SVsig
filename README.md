## *SVsig* - Recurrent structural variations detection in cancer whole genomes


*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/187609)

*MATLAB toolboxes needed:*
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

This version has primarily been tested using `MATLAB_R2020a`. 

##How to run _SVsig_


###Preparing files
SVSig takes in an input .csv file with 15 columns. An example is in this repo, merged_1.6.1.csv. YYour file must match the column names exactly. 
- seqnames, start, strand, altchr, altpos, altstrand: genomic coordinates of both rearrangement breakpoints.
    - Note that chromosome coordinates are integers only. ChrX and ChrY are changed to 23 and 24, respectively. 
- dcc_project_code: histology or tissue type information. 
- sv_id: ID for individual rearrangement
- sid: Patient ID for the rearrangement
- donor_unique_id: Patient ID for the rearrangement
- weights: 

Optional columns:
If this information is not available, set column values to arbitrary value. Will not affect ability to run the model.
- topo: rearrangement topology information
- tppo_n: number of rearrangements involved in topology. 
- mech: DNA damage repair mechanism. 
- homseq: number of base pairs of microhomology at the breakpoint junctions. 
 

###Simple Rearrangements Model (_SVsig-2D_)

SVsig-2D considers each rearrangement to occur independently of each other.
- Open runSVsig.m
  - Change path to sample rearrangements file within lines 22-43
- Open Run2DModel.m 
  - Make sure: complex, weights, and simulation parameters are all false
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run Run2DModel.m

<br>

###Complex Rearrangements Model (_SVsig-2Dc) 
- Open runSVsig.m
  - Change path to sample rearrangements file within lines 22-43
- Open Run2DModel.m
  - Set weights and complex parameters to true. 
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run Run2DModel.m


There are additional parameters that can be set in Run2DModel.m
- 
