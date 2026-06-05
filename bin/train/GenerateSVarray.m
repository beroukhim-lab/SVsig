function [output, Uevent, Usample, Upatient, UTumor, Ustrand1, Ustrand2] = GenerateSVarray(input, length_th, CHR, column_map)
global samp_num

disp('generating events array...');

snowman_table = readtable(input);

if nargin < 4 || isempty(column_map)
    column_map = defaultColumnMap();
end

required_fields = {'seqnames', 'start', 'strand', 'altchr', 'altpos', 'altstrand', 'sid'};
for k = 1:numel(required_fields)
    field_name = required_fields{k};
    if ~isfield(column_map, field_name) || isempty(column_map.(field_name))
        error('GenerateSVarray:MissingColumnMap', ...
              'Missing required column mapping for field: %s', field_name);
    end
end

seqnames_col = get_column_name(snowman_table, column_map.seqnames, true, 'seqnames');
start_col = get_column_name(snowman_table, column_map.start, true, 'start');
strand_col = get_column_name(snowman_table, column_map.strand, true, 'strand');
altchr_col = get_column_name(snowman_table, column_map.altchr, true, 'altchr');
altpos_col = get_column_name(snowman_table, column_map.altpos, true, 'altpos');
altstrand_col = get_column_name(snowman_table, column_map.altstrand, true, 'altstrand');
sid_col = get_column_name(snowman_table, column_map.sid, true, 'sid');

weights_col = get_column_name(snowman_table, column_map.weights, false, 'weights');

n = height(snowman_table);
output = zeros(n, 11);

[Ustrand1, ~, ic_strand1] = unique(snowman_table.(strand_col));
output(:, 3) = ic_strand1;
[Ustrand2, ~, ic_strand2] = unique(snowman_table.(altstrand_col));
output(:, 6) = ic_strand2;

output(:, 2) = to_numeric_column(snowman_table.(start_col), start_col);
output(:, 5) = to_numeric_column(snowman_table.(altpos_col), altpos_col);

UTumor = {'NA'};
output(:, 7) = ones(n, 1);

% Event IDs are internal and not tied to any optional input column.
Uevent = cellstr("event_" + string((1:n)'));
output(:, 8) = (1:n)';

[Usample, ~, ic_sample] = unique(snowman_table.(sid_col));
output(:, 9) = ic_sample;

% Keep sample and patient identity aligned to sid for a single, robust
% identifier path.
Upatient = Usample;
output(:, 10) = output(:, 9);

if ~isempty(weights_col)
    output(:, 11) = to_numeric_column(snowman_table.(weights_col), weights_col);
else
    output(:, 11) = ones(n, 1);
end

% choose samples to include if subsetting and only keep those in output
if ~isempty(samp_num)
    idx = randsample(size(Upatient, 1), samp_num, false);
    idxl = ismember(output(:, 10), idx);
    output = output(idxl, :);
    snowman_table = snowman_table(idxl, :);
end

output(:, 1) = normalize_chr_column(snowman_table.(seqnames_col), seqnames_col);
output(:, 4) = normalize_chr_column(snowman_table.(altchr_col), altchr_col);

% remove lines with events smaller than length_th bp
len_thr_idx = (output(:, 1) == output(:, 4)) & (abs(output(:, 2) - output(:, 5)) < length_th);
output(len_thr_idx, :) = [];

% remove lines with chromosomes out of requested range
for c1 = 1:24
    if sum(CHR == c1) == 0
        chr_idx = (output(:, 1) == c1 | output(:, 4) == c1);
        output(chr_idx, :) = [];
    end
end

end

function column_map = defaultColumnMap()
    column_map = struct( ...
        'seqnames', 'seqnames', ...
        'start', 'start', ...
        'strand', 'strand', ...
        'altchr', 'altchr', ...
        'altpos', 'altpos', ...
        'altstrand', 'altstrand', ...
        'sid', 'sid', ...
        'weights', 'weights');
end

function name = get_column_name(tbl, requested_name, required, label)
    name = '';

    if isempty(requested_name)
        if required
            error('GenerateSVarray:MissingColumn', 'Missing required input column: %s', label);
        end
        return;
    end

    if any(strcmp(tbl.Properties.VariableNames, requested_name))
        name = requested_name;
        return;
    end

    if required
        error('GenerateSVarray:MissingColumn', 'Missing required input column: %s', requested_name);
    end
end

function vals = to_numeric_column(raw_col, col_name)
    if isnumeric(raw_col)
        vals = double(raw_col);
        return;
    end

    if iscell(raw_col)
        raw_col = string(raw_col);
    elseif ischar(raw_col)
        raw_col = string(cellstr(raw_col));
    elseif iscategorical(raw_col)
        raw_col = string(raw_col);
    end

    vals = str2double(string(raw_col));
    if any(isnan(vals))
        error('GenerateSVarray:BadNumericColumn', 'Column %s must contain numeric values.', col_name);
    end
end

function chr_num = normalize_chr_column(raw_col, col_name)
    n_rows = numel(raw_col);
    chr_num = zeros(n_rows, 1);

    for i = 1:n_rows
        chr_num(i) = normalize_chr_value(raw_col(i));
    end

    if any(chr_num < 1 | chr_num > 24 | isnan(chr_num))
        error('GenerateSVarray:BadChromosome', ...
              'Column %s contains noncanonical chromosome values. Use 1-22, X, or Y (with or without chr prefix).', ...
              col_name);
    end
end

function chr_val = normalize_chr_value(raw_value)
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

    text = lower(strtrim(string(raw_value)));
    text = erase(text, "chr");

    if text == "x"
        chr_val = 23;
    elseif text == "y"
        chr_val = 24;
    else
        chr_val = str2double(text);
    end
end
    