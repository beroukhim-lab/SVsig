%%%%%%%%% global variables%%%%%%%%%%%%%%%%%
local = true;
model_exist = false;
len_filter=1e6;
bks_cluster=1;
global FDR_THRESHOLD
FDR_THRESHOLD = 0.1;
global complex
complex = false;                   
global num_breakpoints_per_bin 
num_breakpoints_per_bin=100;
global weights
weights = false;
std_filter = 10;
global simulations
simulations = false; 
global CHR
CHR = 1:23;
global genome_build
genome_build='hg_19';
%hg_19 or hg_38

%%% set working directory%%%%%%
if local
pwd = '/Users/shu/SVsig_labcopy'

else 
pwd = '/xchip_beroukhimlab/Shu/git/SVsig'
end 

global WorkDir
WorkDir = pwd;
addpath(genpath(pwd));

%length_factors=[2, 1.5, 1.25, 1.1, 1.05, 1.025, 1.005];
%length_factors=[1, 1.0001, 1.001, 1.01, 1.05, 1.075, 1.1, 1.2, 1.3, 1.4, 1.5, 2];
%length_factors=[1, 1.0001, 1.001, 1.005, 1.01, 1.025, 1.05, 1.075];
%length_factors=[1, 1.001, 1.01, 1.05, 1.1, 1.15, 1.20, 1.25];
%length_factors=[ 1, 1.25, 1.5, 1.75, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 2];
length_factors=[0.3, 0.31];


%length_factors=[1.12, 1.13, 1.14];


%length_factors=[1, 1.0001];

global length_factor

  
num_pos=zeros(length(length_factors), 2);
for c1 = 1:length(length_factors)
    try
    length_factor=length_factors(c1);
    [hits_table, n_unfiltered, n_filtered] = runSVsig_testlengthdist(length_factors(c1), local, model_exist, len_filter, bks_cluster, FDR_THRESHOLD, complex, num_breakpoints_per_bin, weights, [], std_filter, simulations);
    num_pos(c1, 1)=n_unfiltered;
    num_pos(c1, 2)=n_filtered;
    
    catch ME
    fprintf('no positive results');    
    continue;
    
    end   
end

num_pos1 = [length_factors' num_pos];
label = {'length_factor';'n_unfiltered';'n_filtered'};
pos_results_tble = array2table(num_pos1,'VariableNames',label);


hi=1;


%writetable(pos_results_tble, '/Users/shu/2d_results/20220720_testmarginals_diagsim_withmarg_tilewise_n3_3030bins_Rsum_1e3events_i1','delimiter','\t')

