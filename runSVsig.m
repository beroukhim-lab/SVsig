
function[hits_table] = runSVsig(local, model_exist, len_filter, bks_cluster, FDR_THRESHOLD, complex, num_breakpoints_per_bin, weights, prop_subset, std_filter, simulations)
 
%global samp_num
%samp_num is for the power analysis 
 %%%%%%set data directory%%%%%%%%%%%%%%%%


%DataDir = strcat(WorkDir,'/data/');
%TracksDir = strcat(WorkDir,'/tracks/');

%%%%%%%%% load data table with merged SV with the following columns:
%these are the rearrangments (the events)
% {seqnames, start, strand1, altchr, altpos, strand2,
% subtype(histology)(aka dcc_project_code), sv_id, sid(sample ID),
% donor_unique_id}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ofer's sv_file for ICGC
%sv_file='/Volumes/xchip_beroukhimlab/ofer/matlab/merged_1.6.1.csv';
%my sv_file from Xiatong's adjacency matrix 

if local 
       
    if complex 
 %sv_file ='/Volumes/xchip_beroukhimlab/Kiran/adjancencies/20190502prepped.csv' 
        if weights
%use this first one
%sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/v16prepped_weighted_events_20190724.csv'
 %sv_file='/Users/shu/SVsig_labcopy/complex.csv'
 %sv_file= '/Users/shu/DLBCL/highconf_complex_svaba_1.csv'
 sv_file='/Users/shu/2d_results/v16prepped_weighted_events_20190724.csv'


%sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/v16prepped_weighted_events_20200206.csv'
        end 
       
    elseif simulations 
   
  %sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/20200722simulatedalpha0.csv'
  %sv_file = "/Volumes/xchip_beroukhimlab/Kiran/complex/20200722simulatedalpha1.csv"
  %sv_file = "/Volumes/xchip_beroukhimlab/Kiran/complex/20191202_10nsimulatedalpha1.csv" 
  %sv_file = "/Users/shu/sims/20210204_sij1dx10_ratioslij_a1.csv"
  %sv_file = "/Users/shu/2d_results/2d_results_final/20210527_mixmodel_a0.csv"
  sv_file= '/Users/shu/2d_results/20210722_mixmodel_500kb_mixmodel.csv'
  
    else 
            
 %sv_file='/Volumes/xchip_beroukhimlab/ofer/matlab/merged_1.6.1.csv'
% sv_file='/Volumes/xchip_beroukhimlab/Shu/prepped_intermediate_ccle_subset.csv'
 %sv_file='/Volumes/xchip_beroukhimlab/Shu/merged_1.6.1.csv'
 
 sv_file='/Users/shu/SVsig_labcopy/merged_1.6.1.csv'
 %sv_file='/Users/shu/2d_results/20220803_merged_phgg_combined.csv'
  %sv_file= '/Users/shu/2d_results/20220224_events_greaterthan1e5.csv'
  
  %sv_file='/Users/shu/2d_results/20220315_a1_sim.csv'
  %sv_file='/Users/shu/2d_results/20220316_a1_nolengthdist.csv'
  %sv_file='/Users/shu/2d_results/20220505_newsim_eveninter.csv' %this is
  %the one that gives no results
  

 
    end
    
    
    

else
    if complex
        sv_file = '/xchip/beroukhimlab/Kiran/adjancencies/20190502prepped.csv' 
       if weights 
           sv_file = '/xchip/beroukhimlab/Kiran/complex/v16prepped_weighted_events_20190724.csv'
            sv_file='/Users/shu/DLBCL/highconf_complex_svaba_1.csv'

            %sv_file = '/xchip/beroukhimlab/Kiran/complex/v16prepped_weighted_events_20200206.csv'
       end
        
    elseif simulations 
%%sv_file = "/xchip/beroukhimlab/Kiran/complex/20200722simulatedalpha0.csv"
sv_file = '/xchip/beroukhimlab/Kiran/complex/20200722simulatedalpha1.csv'
%sv_file = "/xchip/beroukhimlab/Kiran/complex/20191202_10nsimulatedalpha1.csv"
%sv_file = "/xchip/beroukhimlab/Kiran/complex/20200422simtoyalpha1.csv"
    else 
 %sv_file ='/xchip/beroukhimlab/ofer/matlab/merged_1.6.1.csv'

    end 
end





SVTable=readtable(sv_file, 'Delimiter', ',');


%%%%%%%%%%load or create background model%%%%%%%%%%%%%%%%%%%

if model_exist
 
    %override bks_cluster assignment in background model loading
    %bks_cluster = 1;
    if complex 
       %load ('backgroundmodel_adjacencies_weighted_20190524');
       load('20210512_complexbackground_5e5bins');
    else 
        %something is not right with this background model?
        %load ('20210421_5e6model_testpval.mat')
        %load 'ICGC_2D_SV_model20190705_88hits.mat'
        %load ('20210423_background_100bkpts.mat')
        %load ('20210423_background_5e5bins.mat')
        %load('20210525_simple.mat')
        load('20210613_mixmodel_a050.mat')
        %load('20210623_mixmodel_a0.mat')
    end 
else
      

    mixmodel;
  %load('/Users/shu/2d_results/20211203_pcawg_background_500kb.mat')
 

   %save('ICGC_2D_SV_model20190705_88hits.mat')
 
end


bks_cluster=0;
if ~bks_cluster
   if weights 
       
       %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom(mfull00, mix_model, [], []);
      % [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom_avgdist(mfull00, mix_model, [], [], bins, events00);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [],1);

    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull+mfull', mix_model, [], [],1, FDR_THRESHOLD);
   else 
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [],1);
   end 

else
    %sij1dx = length_dist_1d_bins(events00,chsize,100);
    %PValMH adjusts for clustered fragile sites within bins whereas PVal does not
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull+mfull', mix_model, bins, events, sij1dx, chsize, CHR, 1, 0.1, 0);
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal_AvgDist(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);

end


%load('/Users/shu/20210916_qfdr_10xpoisson.mat')
%qFDR_mix=qFDR_tophits;
%pa_mix=pa;
%pval_tophits_mix=pval_tophits;
%mfull_pval_mix=mfull;


[hitstable_mix,hitstable_mix_lookup] = HitsTableCV(mfull_pval_mix,pa_mix, pval_tophits_mix, bins_event_tble, qFDR_mix, events00, refgene_tble);

CuratedFusionGene0=CuratedFusionGene(1:end-3,:);
TbyGene_mix = TophitsByGenes(hitstable_mix,hitstable_mix_lookup,1e4,bins,refgene,refgene_tble, [] ,CosmicCencus,CuratedFusionGene0,[]);

%%%%%%%%%FILTRATION%%%%%%%%%%%%%%%%%%%%%%%%%
%only include hits that are interchromosomal or large than the length
%filter

h1=1;
TbyGene_mix_lf = table();
hit_2_include=[];
for c1=1:size(TbyGene_mix,2)
    if (TbyGene_mix(c1).avg_dist == -1 || TbyGene_mix(c1).avg_dist > len_filter)
        hit_2_include(h1)=c1;
        h1=h1+1;
    end
end
TbyGene_mix_lf = TbyGene_mix(hit_2_include);
        


annotated_table = annotate_hits_list( TbyGene_mix_lf,SVTable,bins,hitstable_mix_lookup,pa_mix );
hits_table=table();
hits_table.cluster_num = annotated_table.hit_num;
hits_table.sid = annotated_table.sid;
hits_table.gene_i = annotated_table.gene_i;
hits_table.gene_j = annotated_table.gene_j;

hits_table.subtype = annotated_table.dcc_project_code;
hits_table.chr_i = annotated_table.seqnames;
hits_table.pos_i = annotated_table.start;
hits_table.strand_i = annotated_table.strand;
hits_table.chr_j = annotated_table.altchr;
hits_table.pos_j = annotated_table.altpos;
hits_table.strand_j = annotated_table.altstrand;
hits_table.pval = annotated_table.pval;
hits_table.prob = annotated_table.p_mix;





%annotate number of hits 
%v = accumarray(hits_table.cluster_num, 1);
%for c1= 1:size(hits_table, 1)
%hits_table.num_hits(c1) = v(hits_table.cluster_num(c1));
%end 
    

disp(strcat('the number of hits pre-filtration is ...', num2str(length(unique(hits_table.cluster_num)))))
h1=1;
clusters_to_keep = [];
num_samp = zeros(0,length(unique(hits_table.cluster_num)));
for c1= 1:length(unique(hits_table.cluster_num))

 %find all the samples for the particular cluster_num
 idx = find(hits_table.cluster_num==c1);
  subtable = hits_table(idx, :);
  
  %annotate number of unique samples per hit
 num_samp(c1) = length(unique(subtable.sid));
  
  %are all the hits from the same sample?
  %are the breakpoints greater than a standard deviation of 10?
 if (size(unique(subtable.sid),1) > 1  && std(subtable.pos_i) > std_filter && std(subtable.pos_j) > std_filter)
     cluster_to_keep(h1) = c1; 
     h1 = h1 + 1;

end 
end 

%put num_samp in the hits_table
for c1= 1:size(hits_table, 1)
hits_table.num_hits(c1) = num_samp(hits_table.cluster_num(c1));
end 
    

%keep only the clusters that pass filteration criteria
keep = ismember(hits_table.cluster_num, cluster_to_keep);

hits_table = hits_table(keep, :);

disp(strcat('the number of hits post-filtration is ...', num2str(length(unique(hits_table.cluster_num)))))
%how many SRJs do we find?
%n_hits = length(unique(hits_table.cluster_num));
%writetable(hits_table, '/Volumes/xchip_beroukhimlab/Kiran/adjancencies/power_calculations/20190705hitstableforrobustness.txt','delimiter','\t')

%%%%%%%%%%%FOR POWER ANALYSIS%%%%%%%%%%%%

%create subcurves for power analysis 
% hits_table.prevalent = zeros(size(hits_table, 1),1);
% n_hprev = 0; 
% for c1= 1:length(unique(hits_table.cluster_num))
%    idx = find(hits_table.cluster_num==c1);
%   subtable = hits_table(idx, :);
%   for c2 = 1:length(unique(subtable.subtype))
%   %count the number of unique samples with this particular subtype
%   count = length(unique(subtable.sid));
%   %divide by all samples 
%   total = sum(ismember(SVTable.dcc_project_code,subtable.subtype(c2)));
%   if(count/total > .12) 
%    %disp(num2str(c1))
%    %annotate hitstable
%    hits_table.prevalent(idx) = 1;
%    n_hprev = n_hprev + 1;
%    
%   end 
%   end
% end 


%%%%%%%%%%%%%FOR POWER ANALYSIS%%%%%%%%%%%% 
%read in hits with robustness annotated

% robust_table = readtable('20190708robustness_factor_top_hits.txt');
% [token,remain] = strtok(robust_table.hit_name, "_");
% robust_table.gene_i = token;
% robust_table.gene_j =  erase(remain,"_");
% 
% %merge robust table and hit table on gene_i and gene_j
% robust_hits_table = innerjoin(robust_table, hits_table,'Keys',{'gene_i', 'gene_j'});
% 
% %how many are highly robust?
% idx_high = find(robust_hits_table.p_robust >= 7);
% n_high = length(unique(robust_hits_table(idx_high, :).cluster_num));
% 
% %how many are medium robust?
% idx_med = find(robust_hits_table.p_robust <  7 & robust_hits_table.p_robust >=  2);
% n_med = length(unique(robust_hits_table(idx_med, :).cluster_num));
% 
% 
% %how many are low robust?
% idx_low = find(robust_hits_table.p_robust < 2);
% n_low = length(unique(robust_hits_table(idx_low, :).cluster_num));
end   




