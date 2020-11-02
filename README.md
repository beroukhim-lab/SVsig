## *SVsig* - Recurrent structural variations detection in cancer whole genomes

Description
-----------
*SVsig* is a method developed to classify rearrangements as passenger or driver in cancer patient cohort of whole genome sequences. The distribution of rearrangements in the cancer genome is shaped by both the mechanisms of their formation and the fitness advantages they confer on the cell. This analysis revealed significant predictors of the distribution of rearrangement across the genome and identified known and novel rearrangements that recurred more often than expected given these predictions (for more detailed description: https://doi.org/10.1101/187609)

*MATLAB toolboxes needed:*
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

**To Run _SVsig_ locally **

*Simple Rearrangements Model:*
- Open runSVsig.m
  - Change path to sample rearrangements file within lines 22-43
- Open Run2DModel.m 
  - Keep default parameters
  - Make sure: complex, weights, and simulation parameters are all false
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run Run2DModel.m



*Complex + Simple Rearrangements Model:* 
- Open runSVsig.m
  - Change path to sample rearrangements file within lines 22-43
- Open Run2DModel.m
  - Keep default parameters
  - Set weights and complex parameters to true. 
  - Set working directory in line 48
  - Change path to write hits table in line 83
  - Run Run2DModel.m


