global refgene sij1dx sij1dy

EventsFile=sv_file;

load(strcat(WorkDir,'/tracks/MiscVar')); % run gen_misc_var to update variables
fragile_track=importdata([WorkDir '/tracks/fragile_genes_smith.hg19fp.bed.txt']);

CHR = 1:23; % chromosomes to include in analysis
EventLengthThreshold=200; % filter short event [bp] 
num_breakpoints_per_bin=100; %def= 80
num_annot=4;
min_bin_dist = 500; % def = 500; minimum distance of bins-separting events

Tumor_column=7; 
Event_column=8;
Sample_column=9;
Patient_column=10;

% Generate numeric array of events from merge set
[events0, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2, Utopo, Umech] = GenerateSVarray(EventsFile,EventLengthThreshold,CHR,Tumor_column,Event_column,Sample_column,Patient_column);

% remove events in mask_track
[events0,masked_events] = mask_events( events0,mask_track );
[events0,masked_l1_events] = mask_events( events0,l1_track );

% set bins boundries 
[bins0, numbins] = SetBins(events0,num_breakpoints_per_bin,chsize,CHR,min_bin_dist); % returns a table of bins with chr number, start and end position, and number of breakpoints per bin

% remove bins with low density of events (need to set up threshold manually) 
[bins0, events0, numbins] = remove_low_density_bins(bins0,events0);

[mfull0,bins_event_tble0] = BuildMatrix(events0,bins0,num_annot);

[bins_event_tble, bins, mfull, events, removed_events] = RemoveSameSampleEvents(bins_event_tble0, bins0, mfull0, events0,Patient_column,1);

[bins_event_tble, bins, mfull, events, removed_events_std] = RemoveZeroVarSampleEvents(bins_event_tble, bins, mfull, events);

R = MarginalProbability(bins_event_tble,events,numbins); 

sij1dx = length_dist_1d_bins(events,chsize,10);
sij1dy = EventLengthDist_G(sij1dx,events,0);
sij1dy = sum(sij1dy,2);
sij1dy = sij1dy./sum(sij1dy(1:end-1).*diff(sij1dx'));

annot_array=event_annot(events,TAD_track,fragile_track,gene_track,cancer_genes_track);


no_annot=1;
if no_annot    
    fragile_annot=12;
    sij = ConditionalProbability_fragile(events,annot_array(:,fragile_annot),bins_event_tble,chsize,bins,CHR,sij1dx);
    
    [p, qe, qsolve] = q_solver(R, sij, 1);
    
%    [qFDR_CV, pa_CV, pval_tophits_CV, mfull_pval_CV] = PVal(mfull{1}+mfull{2}+mfull{3}+mfull{4}, p, qe, sij, 1);
else
    [sij,sij1dy] = ConditionalProbability(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx);  % 3D matrix with conditional probability per annotation

    [p, qe, qsolve] = q_solver(R, sum(sij,3), 1);

%    [qFDR_CV, pa_CV, pval_tophits_CV, mfull_pval_CV] = PVal_conv(mfull, qe, sij, 1);
end

%[hitstable_CV,hitstableCV_lookup] = HitsTableCV(mfull_pval_CV,pa_CV, pval_tophits_CV, bins_event_tble, qFDR_CV, events, refgene_tble);

%bins_annot=annotate_bins(bins,bins_event_tble,annot_array,1,l1_track,fragile_track,gene_track,cancer_genes_track,mask_track);

%TbyGene_CV=TophitsByGenes(hitstable_CV,hitstableCV_lookup,1e4,bins,refgene,refgene_tble,UTumor,CosmicCencus,CuratedFusionGene,bins_annot);
