
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
 sv_file ='/Volumes/xchip_beroukhimlab/Kiran/adjancencies/20190502prepped.csv' 
        if weights
sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/v16prepped_weighted_events_20190724.csv'
%sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/v16prepped_weighted_events_20200206.csv'
        end 
       
    elseif simulations 
   
  %sv_file = '/Volumes/xchip_beroukhimlab/Kiran/complex/20200722simulatedalpha0.csv'
  %sv_file = "/Volumes/xchip_beroukhimlab/Kiran/complex/20200722simulatedalpha1.csv"
  %sv_file = "/Volumes/xchip_beroukhimlab/Kiran/complex/20191202_10nsimulatedalpha1.csv" 
  sv_file = "/Users/shu/sims/20210128_binsbydist5e5_noeventratios_sijdx100_a0.csv"
    else 
            
 %sv_file='/Volumes/xchip_beroukhimlab/ofer/matlab/merged_1.6.1.csv'
% sv_file='/Volumes/xchip_beroukhimlab/Shu/prepped_intermediate_ccle_subset.csv'
 %sv_file='/Volumes/xchip_beroukhimlab/Shu/merged_1.6.1.csv'
 sv_file='/Users/shu/SVsig_labcopy/merged_1.6.1.csv'
 %v_file='/Volumes/xchip_beroukhimlab/Shu/ccle_subset_prepped.csv'
 
 
    end
    
    
    

else
    if complex
        sv_file = '/xchip/beroukhimlab/Kiran/adjancencies/20190502prepped.csv' 
       if weights 
           sv_file = '/xchip/beroukhimlab/Kiran/complex/v16prepped_weighted_events_20190724.csv'
            %sv_file = '/xchip/beroukhimlab/Kiran/complex/v16prepped_weighted_events_20200206.csv'
       end
        
    elseif simulations 
%%sv_file = "/xchip/beroukhimlab/Kiran/complex/20200722simulatedalpha0.csv"
sv_file = "/xchip/beroukhimlab/Kiran/complex/20200722simulatedalpha1.csv"
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
       load ('backgroundmodel_adjacencies_weighted_20190524');
    else 
        %something is not right with this background model?
        load 'ICGC_2D_SV_model20190705_88hits.mat'
    end 
  else

    mixmodel;
   %save('ICGC_2D_SV_model20190705_88hits.mat')
    %save('backgroundmodel_adjacencies_weighted_20190524');
   %save('ICGC_2D_SV_model20190610.mat')
end



% EventLengthThreshold=1e2;
%  
% 
% % events array
% %events 0 seems to be an entirely numeric representation of SV table
% events0=zeros(height(SVTable),12);
% 
% if sum(strcmp(SVTable.Properties.VariableNames, 'histology_abbreviation'))>0
%     SVTable.Properties.VariableNames{'histology_abbreviation'} = 'subtype';
% end
% 
% if isa(SVTable.seqnames,'numeric')    
%     events0(:,1)=SVTable.seqnames;
% else
%     events0(:,1)=chr_str2num(SVTable.seqnames)';
% end
% 
% if isa(SVTable.seqnames,'numeric')    
%     events0(:,1)=SVTable.seqnames;
% else
%     events0(:,1)=chr_str2num(SVTable.seqnames)';
% end
% 
% events0(:,2)=SVTable.start;
% [Ustrand1, ia_strand1, ic_strand1]=unique(SVTable.strand);
% events0(:,3)=ic_strand1;
% 
% 
% events0(:,5)=SVTable.altpos;
% [Ustrand2, ia_strand2, ic_strand2]=unique(SVTable.altstrand);
% events0(:,6)=ic_strand2;
% 
% %where is subtype?
% %I'm assuming subtype refers to cancer subtype and is represented by the
% %variable dcc_project_code 
% %[UTumor, ia_code, ic_code]=unique(SVTable.subtype);
% %[UTumor, ia_code, ic_code]=unique(SVTable.dcc_project_code)
% %events0(:,7)=ic_code;
% 
% [Uevent, ia_event, ic_event]=unique(SVTable.sv_id);
% events0(:,8)=ic_event;
% 
% [Usample, ia_sample, ic_sample]=unique(SVTable.sid);
% events0(:,9)=ic_sample;
% 
% if sum(strcmp('donor_unique_id',SVTable.Properties.VariableNames))>0
%     [Upatient, ia_patient, ic_patient]=unique(SVTable.donor_unique_id);
% else
%     [Upatient, ia_patient, ic_patient]=unique(SVTable.sid);
% end
% events0(:,10)=ic_patient;
% 
% %events0(:,11)=SVTable.HOMLEN;
% %events0(:,12)=SVTable.INSLEN;
% 
% disp(strcat('total events from vcfs: ',num2str(length(events0))));
% % filter mask track
% [events0,masked_events] = mask_events( events0,mask_track );
% disp(strcat('total events after masked regions: ',num2str(length(events0))));
% 
% 
%  %events matrix
% %bin the events? %what is happening in this for loop below?
% 
% 
% %%%%%%%%%%%start here if model exists%%%%%%%%%%%%%%%%%%%%%%
% %this matrix is not by rearrangment type (inversion, deletion, duplication)
% %bin the events matrix 
% %using mfull is okay but mfull00 (in break_invasion_model) has more stringent requirements for filtering events 
% mfull0=sparse(length(bins),length(bins));
% bins_event_tble0=zeros(length(events0),3);
% for c1 = 1:length(events0)
%     bini0=find(bins(:,1)==events0(c1,1) & bins(:,2)<=events0(c1,2) & bins(:,3)>=events0(c1,2));
%     binj0=find(bins(:,1)==events0(c1,4) & bins(:,2)<=events0(c1,5) & bins(:,3)>=events0(c1,5));
%     if ~isempty(bini0) && ~isempty(binj0)
%         bins_event_tble0(c1,1) = bini0;
%         bins_event_tble0(c1,2) = binj0;
%         bini = min(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
%         binj = max(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
%         bins_event_tble0(c1,3) = sub2ind([length(bins) length(bins)],bini,binj); % note: assign values only to upper tria of the matrix
%         mfull0(bini,binj) = mfull0(bini,binj) + 1;
%     end
% end
% 
% % remove unassigned events
% unassigned_events=bins_event_tble0(:,1)==0;
% events0(unassigned_events,:)=[];
% bins_event_tble0(unassigned_events,:)=[];
% disp(strcat('total events after removing unassigned events: ',num2str(length(events0))));
% 
% % remove same bin same sample events
% T_bin_sample=table();
% T_bin_sample.bin=bins_event_tble0(:,3);
% T_bin_sample.sample=events0(:,9);
% [u_t,ia_t,ic_t]=unique(T_bin_sample);
% events=events0(ia_t,:);
% bins_event_tble=bins_event_tble0(ia_t,:);
% disp(strcat('total events after removing same sample same bin events: ',num2str(length(events))));
% 
% % remove below length threshold events
% below_length_th=(abs(events(:,5)-events(:,2))<EventLengthThreshold)&(events(:,1)==events(:,4));
% events(below_length_th,:)=[];
% bins_event_tble(below_length_th,:)=[];
% disp(strcat('total events after removing below length threshold: ',num2str(length(events))));
% 
% 
% 
% % update mfull by rebuilding it using events, which is a filtered version
% %of events0 (used to be original mfull0)
% mfull=sparse(length(bins),length(bins));
% for c1 = 1:length(events)
%     bini0=find(bins(:,1)==events(c1,1) & bins(:,2)<=events(c1,2) & bins(:,3)>=events(c1,2));
%     binj0=find(bins(:,1)==events(c1,4) & bins(:,2)<=events(c1,5) & bins(:,3)>=events(c1,5));
%     if ~isempty(bini0) && ~isempty(binj0)
%         mfull(bini0,binj0) = mfull(bini0,binj0) + 1;
%     end  
% end

%in matlab for real numbers mfull' is the transpose of mfull
%use mfull00 over mfull' + mfull because the events in mfull00 are filtered
%more stringently`ewq
if ~bks_cluster
   if weights 
       [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom(mfull00, mix_model, [], []);
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull+mfull', mix_model, [], [],1, FDR_THRESHOLD);
   else 
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [],1);
   end 

else
    sij1dx = length_dist_1d_bins(events00,chsize,10);
    %PValMH adjusts for clustered fragile sites within bins whereas PVal does not
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull+mfull', mix_model, bins, events, sij1dx, chsize, CHR, 1, 0.1, 0);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);
end


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
        

%for export into tile merging, etc. in R 
%dlmwrite('/Volumes/xchip_beroukhimlab/Kiran/adjancencies/20190529wbins.txt', bins, 'delimiter','\t','newline','pc','precision',13);
%dlmwrite('/Volumes/xchip_beroukhimlab/Kiran/adjancencies/20190529whitstable.txt', hitstable_mix, 'delimiter','\t','newline','pc','precision',13);

% %written by Kiran 3/25/19
% %create a list of the unique significant tiles
% sig_tiles = table();
% sig_tiles.chri = bins(hitstable_mix_lookup(:,1), 1); 
% sig_tiles.starti = bins(hitstable_mix_lookup(:,1), 2);
% sig_tiles.endi = bins(hitstable_mix_lookup(:,1), 3);
% sig_tiles.chrj = bins(hitstable_mix_lookup(:,2), 1);
% sig_tiles.startj = bins(hitstable_mix_lookup(:,2), 2);
% sig_tiles.endj = bins(hitstable_mix_lookup(:,2), 3);
% %sig_tile.qval = hitstable_mix(hitstable_mix_lookup(,:1:2),5)

%sig_tiles.numevents = for each row count the number of events that are BOTH in bin i and bin j 
%writetable(unique(sig_tiles), '/xchip/beroukhimlab/Kiran/adjancencies/sigbinstest20190325', 'delimiter','\t')

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
    

%written by Kiran 3/8/19
%remove entries from hits_table that have only supported by 1 sample
% v = accumarray(hits_table.cluster_num, 1);                  % Tally Occurrences Of Rows
% u_cluster_id =  unique(hits_table.cluster_num);  %unique cluster_ids
% indices = u_cluster_id ( v > 1); %Which clusters have greater than 1 occurence
% hits_table = hits_table(ismember(hits_table.cluster_num, indices),:); %kwwp hits have greater than 1 occurence
% disp(strcat(num2str(length(v) - length(v(v>1))), '   hits are supported by only one sample'));

%remove samples that have only 1 unique sample id (i.e. those that are
%supported by only 1 sample)
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




