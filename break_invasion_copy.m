
global refgene sij1dx sij1dy
global num_breakpoints_per_bin
global WorkDir
global complex
global weights
global CHR




EventsFile=sv_file;

load(strcat(WorkDir,'/tracks/MiscVar')); % run gen_misc_var to update variables
fragile_track=importdata([WorkDir '/tracks/fragile_genes_smith.hg19fp.bed.txt']);

% CHR = 1:23; % chromosomes to include in analysis
%now global variable
EventLengthThreshold=200; % filter short event [bp] 
%num_breakpoints_per_bin=1000; %def= 80
num_annot=4;
min_bin_dist = 500; % def = 500; minimum distance of bins-separating events

Tumor_column=7; 
Event_column=8;
Sample_column=9;
Patient_column=10;
Weights_column = 11;

% Generate numeric array of events from merge set
[events0, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2] = GenerateSVarray(EventsFile,EventLengthThreshold,CHR,Tumor_column,Event_column,Sample_column,Patient_column, Weights_column);

%returns an events matrix with a list of junctions
%also returns unique vectors for the sample ids, sv_ids, tumor subtypes,
%strands,topologies and mechanism 

% remove events in mask_track
[events0,masked_events] = mask_events( events0,mask_track );
%returns masked_events which is a logical vector indicating whether or not
%an event is masked?

[events0,masked_l1_events] = mask_events( events0,l1_track );


% set bins boundries 
%[bins0, numbins] = SetBins(events0,num_breakpoints_per_bin,chsize,CHR,min_bin_dist); % returns a table of bins with chr number, start and end position, and number of breakpoints per bin
[bins0, numbins] = SetBins_bylength(events0,5e5, chsize,CHR); % returns a table of bins with chr number, start and end position, and number of breakpoints per bin

% remove bins with low density of events (need to set up threshold manually) 
%UNCOMMENT BELOW -- this throws an error if nothing to remove, fix bug
%later
[bins0, events0, numbins] = remove_low_density_bins(bins0,events0);


%load in some bins
%modelp2=load('/Users/shu/20210121_eventratios_withmarg.mat');
%bins0=modelp2.bins;
%numbins=length(bins0);

[mfull0,bins_event_tble0] = BuildMatrix(events0,bins0,num_annot);

%THIS IS ONLY BC THE NEXT TWO FUNCTIONS ARE COMMENTED OUT
mfull=mfull0;
bins_event_tble=bins_event_tble0;
bins=bins0;
events00=events0;

%UNCOMMENT BELOW
%[bins_event_tble, bins, mfull, events00, removed_events] = RemoveSameSampleEvents(bins_event_tble0, bins0, mfull0, events0,Patient_column,1);

%remove events in the same nucleotide (artifacts)
%UNCOMMENT BELOW
%[bins_event_tble, bins, mfull, events00, removed_events_std] = RemoveZeroVarSampleEvents(bins_event_tble, bins, mfull, events00);

%save bins to compare to TADs size distributions
%dlmwrite(strcat('/Volumes/xchip_beroukhimlab/Kiran/complex/bins', num2str(num_breakpoints_per_bin), '.txt'), bins, 'delimiter','\t','newline','pc','precision',13);

%for export to R to make exploratory graphs
%mfull00 = mfull{1} + mfull{2} + mfull{3} + mfull{4};
%dlmwrite('/Volumes/xchip_beroukhimlab/Kiran/adjancencies/mfull_unweighted_old.txt', nonzeros(mfull00), 'delimiter','\t','newline','pc','precision',13);

R = MarginalProbability(bins_event_tble,events00,numbins); 

%UNCOMMENT 4 ROWS BELOW
sij1dx = length_dist_1d_bins(events00,chsize,1000);
sij1dy = EventLengthDist_G(sij1dx,events00,0);
sij1dy = sum(sij1dy,2);
sij1dy = sij1dy./sum(sij1dy(1:end-1).*diff(sij1dx'));

%sij1dx=repmat(0,100,1);

%no_annot is a boolean indicating whether CP_fragile or CP function should
%be used. What is the difference?
%if complex && weights 
%no_annot = 0; 
%else 
%no_annot=1;
%end 

no_annot = 0; 

if no_annot    
    annot_array=event_annot(events00,TAD_track,fragile_track,gene_track,cancer_genes_track);
    fragile_annot=12;
    sij = ConditionalProbability_fragile(events00,annot_array(:,fragile_annot),bins_event_tble,chsize,bins,CHR,sij1dx);
    
    [p, qe, qsolve] = q_solver_copy(R, sij, 1);
    
%    [qFDR_CV, pa_CV, pval_tophits_CV, mfull_pval_CV] = PVal(mfull{1}+mfull{2}+mfull{3}+mfull{4}, p, qe, sij, 1);
else
    [sij,sij1dy] = ConditionalProbability_copy(events00,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx);  % 3D matrix with conditional probability per annotation

    
    
    mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};
    annot_tiles1=tiles_annot('length',events00,bins,CHR);
     
    short_e=(sum(sum(mfull00(annot_tiles1(:,:,1)))))/(sum(sum(mfull00)));
    long_e=(sum(sum(mfull00(annot_tiles1(:,:,2)))))/(sum(sum(mfull00)));
    inter_e=(sum(sum(mfull00(annot_tiles1(:,:,3)))))/(sum(sum(mfull00)));
    [sij1]=renormalize_tiles(sum(sij,3), short_e, long_e, inter_e, events00, bins, CHR);

    %sij1=sum(sij,3);
    [p, qe, qsolve] = q_solver(R, sum(sij,3), 1);

    short_p=(sum(sum(p(annot_tiles1(:,:,1)))))/(sum(sum(p)))
    long_p=(sum(sum(p(annot_tiles1(:,:,2)))))/(sum(sum(p)))
    inter_p=(sum(sum(p(annot_tiles1(:,:,3)))))/(sum(sum(p)))
    
    %[p]=renormalize_tiles(p, short_e , long_e , inter_e, events00, bins, CHR);

  

%    [qFDR_CV, pa_CV, pval_tophits_CV, mfull_pval_CV] = PVal_conv(mfull, qe, sij, 1);
end

%[hitstable_CV,hitstableCV_lookup] = HitsTableCV(mfull_pval_CV,pa_CV, pval_tophits_CV, bins_event_tble, qFDR_CV, events, refgene_tble);

%bins_annot=annotate_bins(bins,bins_event_tble,annot_array,1,l1_track,fragile_track,gene_track,cancer_genes_track,mask_track);

%TbyGene_CV=TophitsByGenes(hitstable_CV,hitstableCV_lookup,1e4,bins,refgene,refgene_tble,UTumor,CosmicCencus,CuratedFusionGene,bins_annot);
