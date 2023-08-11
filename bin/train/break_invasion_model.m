
global refgene sij1dx sij1dy
global num_breakpoints_per_bin
global WorkDir
global complex
global weights
global CHR
global bins_event_tble
global genome_build
global bin_length



EventsFile=sv_file;

load(strcat(WorkDir,'/data/tracks/MiscVar')); % run gen_misc_var to update variables
fragile_track=importdata([WorkDir '/data/tracks/fragile_genes_smith.hg19fp.bed.txt']);
%fragile_track=importdata([WorkDir 'data//tracks/HumCFS_hg38_lift_hg19.txt']);


if genome_build == 'hg_38'
    fragile_track=importdata([WorkDir '/data/tracks/fragile_genes.hg38.txt']);
    chsize=importdata([WorkDir '/data/tracks/chsize_hg38.txt']);
end

% CHR = 1:23; % chromosomes to include in analysis
%now global variable
EventLengthThreshold=200; % filter short event [bp], deafult 200
%num_breakpoints_per_bin=1000; %def= 80
num_annot=4;
min_bin_dist = 500; % def = 500; minimum distance of bins-separating events

Tumor_column=7; 
Event_column=8;
Sample_column=9;
Patient_column=10;
Weights_column = 11;

% Generate numeric array of events from merge set
%returns an events matrix with a list of junctions
%also returns unique vectors for the sample ids, sv_ids, tumor subtypes, strands,topologies and mechanism 
[events0, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2] = GenerateSVarray(EventsFile,EventLengthThreshold,CHR,Tumor_column,Event_column,Sample_column,Patient_column, Weights_column);

% remove events in mask_track
[events0,masked_events] = mask_events(events0, mask_track);
%returns masked_events which is a logical vector indicating whether or not
%an event is masked?

[events0,masked_l1_events] = mask_events(events0, l1_track);


% set bins boundries 
% returns a table of bins with chr number, start and end position, and number of breakpoints per bin
%[bins0, numbins] = SetBins(events0,num_breakpoints_per_bin,chsize,CHR,min_bin_dist); 
[bins0, numbins] = SetBins_bylength(events0, bin_length, chsize,CHR); 

%bins_all=bins0;

% remove bins with low density of events (need to set up threshold manually) 
%this throws an error if nothing to remo%ve, fix bug later 
[bins0, events0, numbins] = remove_low_density_bins(bins0, events0);

%load in some bins
%modelp2=load('/Users/shu/2d_results/20210722_mixmodel_500kb.mat');
%bins0=modelp2.bins;
%numbins=length(bins0);




%keeping bins with low density but without the eventsi n those bins
%bins0=bins_all;
%numbins=length(bins0);

%bins=table2array(readtable('/Users/shu/2d_results/20220527_1_22_binssub.csv'));
[mfull0,bins_event_tble0] = BuildMatrix(events0, bins0, num_annot);

bins_event_tble=bins_event_tble0;
bins=bins0;
events00=events0;
mfull=mfull0;

[bins_event_tble, bins, mfull, events00, removed_events] = RemoveSameSampleEvents(bins_event_tble0, bins0, mfull0, events0,Patient_column,1);

%remove events in the same nucleotide (artifacts)
[bins_event_tble, bins, mfull, events00, removed_events_std] = RemoveZeroVarSampleEvents(bins_event_tble, bins, mfull, events00);



R = MarginalProbability(bins_event_tble,events00,numbins); 


% 
sij1dx = length_dist_1d_bins(events00,chsize,10);
sij1dx=unique(sij1dx);
sij1dy = EventLengthDist_G(sij1dx,events00,0);
sij1dy = sum(sij1dy,2);
sij1dy = sij1dy./sum(sij1dy(1:end-1).*diff(sij1dx'));

%sij1dx=repmat(0,100,1);

%no_annot indicates if CP_fragile or CP is to be used 
%default is CP
no_annot = 0; 

if no_annot    
    annot_array=event_annot(events00,TAD_track,fragile_track,gene_track,cancer_genes_track);
    fragile_annot=12;
    sij = ConditionalProbability_fragile(events00,annot_array(:,fragile_annot),bins_event_tble,chsize,bins,CHR,sij1dx);
    
    [p, qe, qsolve] = q_solver_copy(R, sij, 1);
    
else
    [sij,sij1dy] = ConditionalProbability_copy(events00,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx);  % 3D matrix with conditional probability per annotation

    
    mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};
    annot_tiles1=tiles_annot_copy('length',events00,bins,CHR);
     
    %normalize by event ratios
    [sij1]=renormalize_tiles(mfull00, sum(sij,3), events00, bins, CHR);

    [p, qe, qsolve] = q_solver(R, sij1, 1);
    
    %normalize again by event ratios 
    [p]=renormalize_tiles(mfull00, p,  events00, bins, CHR);

    %make sure sum to 2
    p = 2*p ./ sum(sum(p));


end

