
global refgene sij1dx sij1dy
global num_breakpoints_per_bin
global WorkDir
global complex
global weights
global CHR
global bins_event_tble
global genome_build
global bin_length
global BIN_SHIFT

global TRACK_PATHS
global LOW_DENSITY_THRESHOLD
global refgene_lookup
global cancer_gene_symbols
global curated_fusion_pairs
global chromosome_sizes

EventsFile=sv_file;

% Load track files using resolved paths from run2DModel
tracks = load_track_files(TRACK_PATHS, genome_build);
refgene = tracks.refgene;
refgene_lookup = tracks.refgene_lookup;
chromosome_sizes = tracks.chromosome_sizes;
cancer_gene_symbols = tracks.cancer_gene_symbols;
curated_fusion_pairs = tracks.curated_fusion_pairs;
blacklist_regions = tracks.blacklist_regions;
l1_elements = tracks.l1_elements;

% CHR = 1:23; % chromosomes to include in analysis
%now global variable
EventLengthThreshold=200; % filter short event [bp], deafult 200
%num_breakpoints_per_bin=1000; %def= 80
num_annot=4;
min_bin_dist = 500; % def = 500; minimum distance of bins-separating events

sample_id_column = 9;

input_columns = struct( ...
	'seqnames', 'seqnames', ...
	'start', 'start', ...
	'strand', 'strand', ...
	'altchr', 'altchr', ...
	'altpos', 'altpos', ...
	'altstrand', 'altstrand', ...
	'sid', 'sid', ...
	'weights', 'weights');

% Generate numeric array of events from merge set
%returns an events matrix with a list of junctions
%also returns unique vectors for the sample ids, sv_ids, tumor subtypes, strands,topologies and mechanism 
[events0, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2] = GenerateSVarray(EventsFile, EventLengthThreshold, CHR, input_columns);

% remove events in blacklist regions
[events0,masked_events] = mask_events(events0, blacklist_regions);
%returns masked_events which is a logical vector indicating whether or not
%an event is masked?

[events0,masked_l1_events] = mask_events(events0, l1_elements);


% set bins boundries 
% returns a table of bins with chr number, start and end position, and number of breakpoints per bin
%[bins0, numbins] = SetBins(events0,num_breakpoints_per_bin,chromosome_sizes,CHR,min_bin_dist);
[bins0, numbins] = SetBins_bylength(events0, bin_length, chromosome_sizes, CHR, BIN_SHIFT);

% 2026/02/17: exported bins0 as the bin boundary definitions, before we
% drop the empty bins

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

[bins_event_tble, bins, mfull, events00, removed_events] = RemoveSameSampleEvents(bins_event_tble0, bins0, mfull0, events0, sample_id_column, 1);

%remove events in the same nucleotide (artifacts)
[bins_event_tble, bins, mfull, events00, removed_events_std] = RemoveZeroVarSampleEvents(bins_event_tble, bins, mfull, events00);

% Keep the user-requested chromosome scope and report which of those were dropped.
requested_CHR = CHR;
% Antonia: Update CHR to only include bins with at least one event
CHR = unique(bins(:,1))';
excluded = setdiff(requested_CHR, CHR);
disp(['Excluded chromosomes: ' num2str(excluded)])

R = MarginalProbability(bins_event_tble,events00,numbins); 


% 
sij1dx = length_dist_1d_bins(events00,chromosome_sizes,10);
sij1dx=unique(sij1dx);
sij1dy = EventLengthDist_G(sij1dx,events00,0);
sij1dy = sum(sij1dy,2);
sij1dy = sij1dy./sum(sij1dy(1:end-1).*diff(sij1dx'));

%sij1dx=repmat(0,100,1);

% remove CP fragile as an option
[sij,sij1dy] = ConditionalProbability(events00,chromosome_sizes,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx);  % 3D matrix with conditional probability per annotation


mfull00=mfull{1}+mfull{2}+mfull{3}+mfull{4};
annot_tiles1=tiles_annot('length',events00,bins,CHR);
    
%normalize by event ratios
[sij1]=renormalize_tiles(mfull00, sum(sij,3), events00, bins, CHR);

[p, qe, qsolve] = q_solver(R, sij1, 1);

%normalize again by event ratios 
[p]=renormalize_tiles(mfull00, p,  events00, bins, CHR);

%make sure sum to 2
p = 2*p ./ sum(sum(p));

