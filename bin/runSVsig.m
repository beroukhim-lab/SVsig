
function[hits_table] = runSVsig(sv_file, model_exist, complex, weights, len_filter, bks_cluster, ...
FDR_THRESHOLD, bin_length, num_breakpoints_per_bin, std_filter)
 

%load rearrangement data table with the following columns:

% {seqnames, start, strand1, altchr, altpos, strand2,
% subtype(histology), sv_id, sid(sample ID),donor_unique_id}

SVTable=readtable(sv_file, 'Delimiter', ',');



%%%%%%%%%%load or create background model%%%%%%%%%%

if model_exist
 
    %override bks_cluster assignment in background model loading
    %bks_cluster = 1;
    if complex 
       %load ('backgroundmodel_adjacencies_weighted_20190524');
       load('20210512_complexbackground_5e5bins');
    else 
        %load('20210613_mixmodel_a050.mat')
        load('20240719_testmodel.mat')
    end 
else
      
    %calculate break invasion and double break join models
    %get optimal combination of both models: mix_model 
    mixmodel;

    
end


bks_cluster=1;
if ~bks_cluster
   if weights 
       %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom(mfull00, mix_model, [], []);
       %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom_avgdist(mfull00, mix_model, [], [], bins, events00);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [],1);

   else 
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull00, mix_model, [], [],1);
   end 

else
    sij1dx = length_dist_1d_bins(events00,chsize,100);
    %PValMH adjusts for clustered fragile sites within bins whereas PVal does not
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull+mfull', mix_model, bins, events, sij1dx, chsize, CHR, 1, 0.1, 0);
    %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValMH(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal_AvgDist(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);
     %[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal_AvgDist_copy(mfull00, mix_model, bins, events00, sij1dx, chsize, CHR, 1, 0);

end




%this function can take a long time bc loops through all possible genes and
%hits and such 
%wait jk TophistByGenes doe sthat check HitsTableCV later 

[hitstable_mix,hitstable_mix_lookup] = HitsTableCV(mfull_pval_mix,pa_mix, pval_tophits_mix, bins_event_tble, qFDR_mix, events00, refgene_tble);

disp ('done with HitsTableCV')

CuratedFusionGene0=CuratedFusionGene(1:end-3,:);
TbyGene_mix = TophitsByGenes(hitstable_mix,hitstable_mix_lookup,1e4,bins,refgene,refgene_tble, [] ,CosmicCencus,CuratedFusionGene0,[]);




%%%%%%%%%FILTRATION%%%%%%%%%%
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
    

%keep only the clusters that pass filtration criteria
keep = ismember(hits_table.cluster_num, cluster_to_keep);
hits_table = hits_table(keep, :);

disp(strcat('the number of hits post-filtration is ...', num2str(length(unique(hits_table.cluster_num)))))

end   




