%%%Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu
%%% Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
%%% Date last updated: August 11, 2023 

%%% License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute
%%% Dependencies: 
%%% See README for guide on how to run this package

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set working directory %%%
pwd = '/Users/smehrete/Downloads/SVsig/'

global WorkDir
WorkDir = pwd;
addpath(genpath(pwd));

%path to rearrangements_file
%sv_file='/data/merged_1.6.1.csv'
sv_file='/Users/smehrete/Downloads/SVsig/data/2D_input_CCLE_0627.csv'




%path to write output file with significant 2D hits %
output_file = '/Users/smehrete/SVsig_labcopy/results/07022025_CCLE_tophits_2D_5e5.txt'


%%%%%%%%% set additional parameters %%%%%%%%%

%to use loaded background model set to true
%otherwise false to recalculate background model 
model_exist = false;

%run simple or complex model 
global complex
complex = false;   

%are individual rearrangements weighted. 
%set true if using complex model 
global weights
weights = false;



%@param after running the model we filter out hits that have a standard
%deviation of the following in either the first or second break
std_filter = 10;
%filter out all events that are less than 1MB at the end
len_filter=1e6;
%@param bks_cluster determines whether PVal(fragile sites not accounted for) or PValMH(fragile sites accounted for) is used.
%0 => PVal, 1=> PValMH. PValMH = 1 to replicate ofer's final results for
%the ICGC paper
%get very different results for when PVal = 0
bks_cluster=1;
%param FDR_threshold sets the q value cut off, default 0.1
global FDR_THRESHOLD
FDR_THRESHOLD = 0.1;
%set random seed for reproducibility%
%rng(3014)


global num_breakpoints_per_bin %to determine custom bin boundaries for adjacency matrix
num_breakpoints_per_bin=100;

global bin_length %length of bins for adjacency matrix 
bin_length=5e5;

global CHR %which chromosomes to consider in building the matrix
CHR = 1:23;
global genome_build %hg_19 or hg_38
genome_build='hg_38';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%run model
hits_table = runSVsig(sv_file, model_exist, complex, weights, len_filter, bks_cluster, ...
                      FDR_THRESHOLD, bin_length, num_breakpoints_per_bin, ...
                      std_filter);
%write output 
writetable(hits_table, output_file,'delimiter','\t')


