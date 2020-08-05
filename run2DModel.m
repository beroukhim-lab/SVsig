%%%%%%%%% global variables%%%%%%%%%%%%%%%%%
%@param local  = false means run on local machine, local = true means runs on server 
local = true;
%@param model_exist: to re-run background model set to false, to use loaded background model set to true 
model_exist = false;
%filter out all events that are less than 1MB at the end
len_filter=1e6;
%@param bks_cluster determines whether PVal(fragile sites not accounted for) or PValMH(fragile sites accounted for) is used.
%0 => PVal, 1=> PValMH. PValMH = 1 to replicate ofer's final results for
%the ICGC paper
%get very different results for when PVal = 0
 bks_cluster=1;
 %param FDR_threshold sets the q value cut off
  global FDR_THRESHOLD
 FDR_THRESHOLD = 0.1;
 %set random seed for reproducibility%
%rng(3014)
%@param output_file: sets the name for the output file containing the significant 2D hits %
%output_file = 'hitsalljunctions.txt'
%@param num_breakpoints_per bin determines how many breakpoints are in one
%bin (the smaller this is, the more bins are created and vice versa).
%Adjust based on the size of the dataset! 
%@param complex determines whether to use Ofer's v1.6 simple events or
%Xiatong's jabba complex events 
global complex
complex = false;                   
global num_breakpoints_per_bin %Ofer's default was 100%Kiran's default for Xiaotong's Feb matrix was 1000
%else %300 for complex weighted model
num_breakpoints_per_bin=100;
%end
%@param weights
%Do the events have non-integer weights?
global weights
weights = false;
%@param after running the model we filter out hits that have a standard
%deviation of the following in either the first or second break
std_filter = 10;
%@param simulations, run the model on simulated data (yes == 1, no == 0)
global simulations
simulations = false; 
%which chromosomes to consider in building the matrix
global CHR
CHR = 1:23;

%%% set working directory%%%%%%
if local
%pwd = '/Volumes/xchip_beroukhimlab/Kiran/git/2dmodel/SVsig' 
pwd = '/Volumes/xchip_beroukhimlab/Shu/git/SVsig'
else 
%pwd = '/xchip/beroukhimlab/Kiran/git/2dmodel/SVsig'
pwd = '/xchip_beroukhimlab/Shu/git/SVsig'
end 

global WorkDir
WorkDir = pwd;
addpath(genpath(pwd));


%run at different len_filters and number of breakpoints per bin at double the background rate
%num_hits = zeros(2, 4);
%for c1 = 1:2
 %for c2 = 1:4
%lf = [0, 1e6]; 
%nbb = [100, 200, 500, 1000]; 

%run model
%num_breakpoints_per_bin = nbb(c2);
%disp(num_breakpoints_per_bin)
%hits_table = runSVsig(local, model_exist, lf(c1), bks_cluster, FDR_THRESHOLD, complex, num_breakpoints_per_bin, weights, [], std_filter, simulations);
%num_hits(c1, c2) = length(unique(hits_table.cluster_num));
%writetable(hits_table, strcat('/Volumes/xchip_beroukhimlab/Kiran/complex/20191113simulationsbin', num2str(num_breakpoints_per_bin),'filter', num2str(lf(c1)), '.txt'),'delimiter','\t')

 
%end 
%end 

%run model

%note to self: need to check complex, weights, simulations; make sure path
%to appropriate file is correct in "runSVsig.m"
hits_table = runSVsig(local, model_exist, len_filter, bks_cluster, FDR_THRESHOLD, complex, num_breakpoints_per_bin, weights, [], std_filter, simulations);
%writetable(hits_table, '/xchip/beroukhimlab/Kiran/complex/20200212regularbackgroundrate5e5lengthfilter','delimiter','\t')
writetable(hits_table, '/Volumes/xchip_beroukhimlab/Shu/ccle2_sample_practice','delimiter','\t')

%note to self: commented out tile_density line 174 to see what happens 

