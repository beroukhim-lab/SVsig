
function [hits_table, bins] = runSVsig(sv_file, model_exist, complex, weights, len_filter, ...
FDR_THRESHOLD, bin_length, num_breakpoints_per_bin, std_filter, model_file, tier_std_cutoff)
global cancer_gene_symbols
global curated_fusion_pairs
global chromosome_sizes
global refgene
global refgene_lookup
global events00
global bins
global bins_event_tble
global mfull00
global mix_model

if nargin < 10 || isempty(model_file)
    model_file = '';
end
if nargin < 11 || isempty(tier_std_cutoff)
    tier_std_cutoff = 38491;
end
 

% Load rearrangement data table by required column names.
% Required: seqnames, start, strand, altchr, altpos, altstrand,
%           sid

SVTable=readtable(sv_file, 'Delimiter', ',');
SVTable = standardizeSVInputTable(SVTable);



%%%%%%%%%%load or create background model%%%%%%%%%%

if model_exist

    if isempty(model_file)
        error('runSVsig:MissingModelFile', 'model_exist=true requires a non-empty model_file path.');
    end
    if ~isfile(model_file)
        error('runSVsig:MissingModelFile', 'model_file does not exist: %s', model_file);
    end
    load(model_file);
else
      
    %calculate break invasion and double break join models
    %get optimal combination of both models: mix_model 
    mixmodel;

    
end


if complex
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PValCBinom(mfull00, mix_model, [], []);
else
    sij1dx = length_dist_1d_bins(events00,chromosome_sizes,100);
    [qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal_AvgDist(mfull00, mix_model, bins, events00, 0);
end

if isempty(pval_tophits_mix) || size(pval_tophits_mix,1) == 0
    % if no significant hits, return empty table with correct columns and skip annotation steps
    fprintf('[runSVsig] No significant hits found. Writing empty output table.\n');
    hits_table = createEmptyHitsTable();
    return;
end




%this function can take a long time bc loops through all possible genes and
%hits and such 
%wait jk TophistByGenes doe sthat check HitsTableCV later 

[hitstable_mix,hitstable_mix_lookup] = HitsTableCV(mfull_pval_mix,pa_mix, pval_tophits_mix, bins_event_tble, qFDR_mix, events00, refgene_lookup);

disp ('done with HitsTableCV')

TbyGene_mix = TophitsByGenes(hitstable_mix,hitstable_mix_lookup,1e4,bins,refgene,refgene_lookup, [] ,cancer_gene_symbols,curated_fusion_pairs,[]);




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
        


annotated_table = annotate_hits_list( TbyGene_mix_lf,SVTable,bins,hitstable_mix_lookup,pa_mix,qFDR_mix );
hits_table=table();
hits_table.cluster_num = annotated_table.hit_num;
hits_table.sid = annotated_table.sid;
hits_table.gene_i = annotated_table.gene_i;
hits_table.all_genes_i = annotated_table.nearby_genes_i;
hits_table.gene_j = annotated_table.gene_j;
hits_table.all_genes_j = annotated_table.nearby_genes_j;

hits_table.bin_i = annotated_table.bin_i;
hits_table.chr_i = annotated_table.seqnames;
hits_table.pos_i = annotated_table.start;
hits_table.strand_i = annotated_table.strand;
hits_table.bin_j = annotated_table.bin_j;
hits_table.chr_j = annotated_table.altchr;
hits_table.pos_j = annotated_table.altpos;
hits_table.strand_j = annotated_table.altstrand;
hits_table.sv_class = classifySVClass(hits_table.chr_i, hits_table.chr_j, hits_table.strand_i, hits_table.strand_j);
hits_table.tile_qval = annotated_table.tile_qval;
hits_table.pval = annotated_table.pval;
hits_table.prob = annotated_table.p_mix;




disp(strcat('the number of hits pre-filtration is ...', num2str(length(unique(hits_table.cluster_num)))))
clusters = unique(hits_table.cluster_num);
num_clusters = numel(clusters);
clusters_to_keep = [];
num_samp = zeros(num_clusters,1);
stddev_i_by_cluster = zeros(num_clusters,1);
stddev_j_by_cluster = zeros(num_clusters,1);

for c1 = 1:num_clusters

% find all the samples for the particular cluster_num
 idx = find(hits_table.cluster_num == clusters(c1));
 subtable = hits_table(idx, :);

 % annotate number of unique samples per hit
 num_samp(c1) = length(unique(subtable.sid));

 % calculate stddev for both breakpoints
 stddev_i_by_cluster(c1) = std(subtable.pos_i);
 stddev_j_by_cluster(c1) = std(subtable.pos_j);

    % are all the hits from the same sample?
    % are the breakpoints greater than a standard deviation of 10?
 if (size(unique(subtable.sid),1) > 1 && stddev_i_by_cluster(c1) > std_filter && stddev_j_by_cluster(c1) > std_filter)
     clusters_to_keep(end+1) = clusters(c1);
 end
end

[~, cluster_idx] = ismember(hits_table.cluster_num, clusters);
hits_table.num_hits = num_samp(cluster_idx);
hits_table.stddev_i = stddev_i_by_cluster(cluster_idx);
hits_table.stddev_j = stddev_j_by_cluster(cluster_idx);

tier = repmat({'tier 3'}, height(hits_table), 1);
is_tier1 = hits_table.stddev_i <= tier_std_cutoff & hits_table.stddev_j <= tier_std_cutoff;
is_tier2 = xor(hits_table.stddev_i <= tier_std_cutoff, hits_table.stddev_j <= tier_std_cutoff);
tier(is_tier1) = {'tier 1'};
tier(is_tier2) = {'tier 2'};
hits_table.tier = tier;
    

%keep only the clusters that pass filtration criteria
keep = ismember(hits_table.cluster_num, clusters_to_keep);
hits_table = hits_table(keep, :);

disp(strcat('the number of hits post-filtration is ...', num2str(length(unique(hits_table.cluster_num)))))

end   

function SVTable = standardizeSVInputTable(SVTable)
required_cols = {'seqnames', 'start', 'strand', 'altchr', 'altpos', 'altstrand', 'sid'};
for k = 1:numel(required_cols)
    if ~any(strcmp(SVTable.Properties.VariableNames, required_cols{k}))
        error('runSVsig:MissingInputColumn', 'Missing required input column: %s', required_cols{k});
    end
end

SVTable.seqnames = normalizeChrVector(SVTable.seqnames, 'seqnames');
SVTable.altchr = normalizeChrVector(SVTable.altchr, 'altchr');
SVTable.start = toNumericVector(SVTable.start, 'start');
SVTable.altpos = toNumericVector(SVTable.altpos, 'altpos');
end

function values = toNumericVector(raw_values, column_name)
if isnumeric(raw_values)
    values = double(raw_values);
    return;
end

if iscategorical(raw_values)
    raw_values = string(raw_values);
elseif iscell(raw_values)
    raw_values = string(raw_values);
end

values = str2double(string(raw_values));
if any(isnan(values))
    error('runSVsig:BadNumericInput', 'Column %s must contain numeric values.', column_name);
end
end

function chr_num = normalizeChrVector(raw_values, column_name)
n_rows = numel(raw_values);
chr_num = zeros(n_rows, 1);
for i = 1:n_rows
    chr_num(i) = normalizeChrValue(raw_values(i));
end

if any(chr_num < 1 | chr_num > 24 | isnan(chr_num))
    error('runSVsig:BadChromosomeInput', ...
          'Column %s contains noncanonical chromosome values. Use 1-22, X, or Y (with or without chr prefix).', ...
          column_name);
end
end

function chr_val = normalizeChrValue(raw_value)
if iscell(raw_value)
    raw_value = raw_value{1};
end

if isnumeric(raw_value)
    chr_val = double(raw_value);
    return;
end

if iscategorical(raw_value)
    raw_value = string(raw_value);
end

chr_text = lower(strtrim(string(raw_value)));
chr_text = erase(chr_text, "chr");

if chr_text == "x"
    chr_val = 23;
elseif chr_text == "y"
    chr_val = 24;
else
    chr_val = str2double(chr_text);
end
end

function hits_table = createEmptyHitsTable()
var_names = {'cluster_num','sid','gene_i','all_genes_i','gene_j','all_genes_j', ...
             'bin_i','chr_i','pos_i','strand_i','bin_j','chr_j','pos_j','strand_j', ...
             'sv_class','tile_qval','pval','prob','num_hits','stddev_i','stddev_j','tier'};
var_types = {'double','cell','cell','cell','cell','cell', ...
             'double','double','double','double','double','double','double','double', ...
             'cell','double','double','double','double','double','double','cell'};

hits_table = table('Size',[0 numel(var_names)], 'VariableTypes', var_types, 'VariableNames', var_names);
end

function sv_class = classifySVClass(chr_i, chr_j, strand_i, strand_j)
n_rows = numel(chr_i);
sv_class = repmat({''}, n_rows, 1);

for i = 1:n_rows
    if chr_i(i) ~= chr_j(i)
        sv_class{i} = 'inter_chr';
        continue;
    end

    si = normalizeStrandValue(strand_i(i));
    sj = normalizeStrandValue(strand_j(i));

    if si == "unknown" || sj == "unknown"
        sv_class{i} = 'unknown';
    elseif si == "+" && sj == "-"
        sv_class{i} = 'del';
    elseif si == "-" && sj == "+"
        sv_class{i} = 'tandem_dup';
    elseif (si == "+" && sj == "+") || (si == "-" && sj == "-")
        sv_class{i} = 'inv';
    else
        sv_class{i} = 'unknown';
    end
end
end

function strand_text = normalizeStrandValue(raw_strand)
if iscell(raw_strand)
    raw_strand = raw_strand{1};
end

if isnumeric(raw_strand)
    if isscalar(raw_strand)
        if raw_strand == 1
            strand_text = "+";
            return;
        elseif raw_strand == -1
            strand_text = "-";
            return;
        end
    end
    strand_text = "unknown";
    return;
end

if iscategorical(raw_strand)
    raw_strand = string(raw_strand);
end

strand_text = string(raw_strand);
strand_text = strtrim(strand_text);

if strand_text == "+" || strand_text == "-"
    return;
end

if strand_text == "1" || lower(strand_text) == "+1"
    strand_text = "+";
    return;
end

if strand_text == "-1"
    strand_text = "-";
    return;
end

strand_text = "unknown";
end




