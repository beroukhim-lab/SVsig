function tracks = load_track_files(resolved_track_paths, genome_build)
    % Load and parse genomic track files for SVsig analysis.
    %
    % Args:
    %   resolved_track_paths: struct with fields (from resolve_track_paths):
    %       ref_genes_file, cancer_genes_file, curated_fusions_file,
    %       chrom_sizes_file, blacklist_file, l1_elements_file
    %   genome_build: 'hg_19' or 'hg_38'
    %
    % Returns:
    %   tracks: struct with parsed data:
    %       refgene (struct), refgene_lookup (matrix), 
    %       cancer_gene_symbols (cell), curated_fusion_pairs (cell),
    %       chromosome_sizes (vector), blacklist_regions (matrix), l1_elements (matrix)
    
    fprintf('[load_track_files] Loading track files for %s\n', genome_build);
    
    % Load chromosome sizes first (needed for other parsers)
    tracks.chromosome_sizes = load_chrom_sizes(resolved_track_paths.chrom_sizes_file, genome_build);
    fprintf('[load_track_files] Loaded chrom.sizes: %d chromosomes\n', numel(tracks.chromosome_sizes));
    
    % Load reference genes (creates both refgene struct and refgene_lookup)
    [tracks.refgene, tracks.refgene_lookup] = load_ref_genes(resolved_track_paths.ref_genes_file, tracks.chromosome_sizes, genome_build);
    fprintf('[load_track_files] Loaded reference genes: %d entries\n', length(tracks.refgene));
    
    % Load cancer genes (simple list)
    tracks.cancer_gene_symbols = load_cancer_genes(resolved_track_paths.cancer_genes_file);
    fprintf('[load_track_files] Loaded cancer genes: %d entries\n', length(tracks.cancer_gene_symbols));
    
    % Load curated fusions (pairs)
    tracks.curated_fusion_pairs = load_curated_fusions(resolved_track_paths.curated_fusions_file);
    fprintf('[load_track_files] Loaded curated fusions: %d pairs\n', size(tracks.curated_fusion_pairs, 1));
    
    % Load blacklist regions
    tracks.blacklist_regions = load_bed_file(resolved_track_paths.blacklist_file, tracks.chromosome_sizes, genome_build);
    fprintf('[load_track_files] Loaded blacklist regions: %d entries\n', size(tracks.blacklist_regions, 1));
    
    % Load L1 elements
    tracks.l1_elements = load_bed_file(resolved_track_paths.l1_elements_file, tracks.chromosome_sizes, genome_build);
    fprintf('[load_track_files] Loaded L1 elements: %d entries\n', size(tracks.l1_elements, 1));
    
end

function chsize = load_chrom_sizes(chrom_sizes_file, genome_build)
    % Load chromosome sizes and map to numeric vector [1..24].
    % Preferred columns are chr and size; headerless two-column files are also accepted.
    
    validate_file_for_genome_build(chrom_sizes_file, genome_build);

    data = read_chr_size_table(chrom_sizes_file);
    
    chsize = zeros(1, 24);  % chr1 to chrX/Y (24 canonical chromosomes)
    dropped_rows = 0;
    
    for i = 1:height(data)
        chr_name = data.chr{i};
        size_val = data.size(i);
        
        chr_num = parse_chr_name(chr_name);
        if is_canonical_chr_num(chr_num)
            chsize(chr_num) = size_val;
        else
            dropped_rows = dropped_rows + 1;
        end
    end

    if dropped_rows > 0
        fprintf('[load_track_files] Ignored %d noncanonical rows in chrom_sizes_file\n', dropped_rows);
    end
end

function [refgene, refgene_lookup] = load_ref_genes(ref_genes_file, chromosome_sizes, genome_build)
    % Load reference genes using canonical genomic coordinates.
    % Required logical columns are:
    %   chr: chromosome name (aliases: chrom, chromosome)
    %   start: gene start (aliases: txStart, chromStart)
    %   end: gene end (aliases: txEnd, chromEnd, xEnd)
    %   gene_symbol: symbol (aliases: name2, gene, symbol)
    
    validate_file_for_genome_build(ref_genes_file, genome_build);

    data = readtable(ref_genes_file, 'FileType', 'text', 'Delimiter', '\t', ...
                    'ReadVariableNames', true);

    chr_col = get_required_column(data, {'chr', 'chrom', 'chromosome'}, ref_genes_file, 'ref_genes_file');
    start_col = get_required_column(data, {'start', 'txStart', 'chromStart'}, ref_genes_file, 'ref_genes_file');
    end_col = get_required_column(data, {'end', 'txEnd', 'chromEnd', 'xEnd'}, ref_genes_file, 'ref_genes_file');
    gene_symbol_col = get_required_column(data, {'gene_symbol', 'name2', 'gene', 'symbol'}, ref_genes_file, 'ref_genes_file');
    
    n = height(data);
    refgene = struct('rg', struct('locus_id', {}, 'symb', {}, 'start', {}, 'end', {}));
    rg_list = [];
    refgene_lookup = [];
    locus_idx = 0;
    dropped_rows = 0;
    
    for i = 1:n
        chr_name = data{i, chr_col}{1};
        chr_num = parse_chr_name(chr_name);
        start_pos = data{i, start_col};
        end_pos = data{i, end_col};
        symbol = data{i, gene_symbol_col}{1};
        
        if is_canonical_chr_num(chr_num)
            locus_idx = locus_idx + 1;
            rg_list = [rg_list; struct('locus_id', locus_idx, 'symb', symbol, ...
                                       'start', start_pos, 'end', end_pos)];
            % Add to lookup table: [chr_num, start, end, locus_id]
            refgene_lookup = [refgene_lookup; chr_num, start_pos, end_pos, locus_idx];
        else
            dropped_rows = dropped_rows + 1;
        end
    end
    
    refgene.rg = rg_list;

    if dropped_rows > 0
        fprintf('[load_track_files] Ignored %d noncanonical rows in ref_genes_file\n', dropped_rows);
    end
end

function cancer_genes = load_cancer_genes(cancer_genes_file)
    % Load cancer gene symbols.
    % Expected canonical column is gene_symbol; if absent, the first column is used.
    
    try
        data = readtable(cancer_genes_file, 'FileType', 'text', 'Delimiter', {'\t', ','});
        gene_symbol_col = find_column_by_alias(data, {'gene_symbol', 'gene', 'symbol', 'Gene Symbol'});
        if isempty(gene_symbol_col)
            gene_symbol_col = 1;
        end
        cancer_genes = table2cell(data(:, gene_symbol_col));
    catch
        % Fallback: read as plain text, one gene per line
        fid = fopen(cancer_genes_file, 'r');
        cancer_genes = {};
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && ~isempty(line) && ~startsWith(line, '#')
                cancer_genes{end+1} = strtrim(line);
            end
        end
        fclose(fid);
    end
end

function curated_fusions = load_curated_fusions(curated_fusions_file)
    % Load curated fusion pairs
    % Expected canonical columns are gene_A and gene_B.
    % Accepted aliases: gene1/5prime_gene and gene2/3prime_gene.
    
    data = readtable(curated_fusions_file, 'FileType', 'text', 'Delimiter', '\t', ...
                    'ReadVariableNames', true);
    
    % Find gene A column (try multiple names)
    gene_a_col = get_required_column(data, {'gene_A', 'gene1', '5prime_gene'}, curated_fusions_file, 'curated_fusions_file');
    gene_b_col = get_required_column(data, {'gene_B', 'gene2', '3prime_gene'}, curated_fusions_file, 'curated_fusions_file');
    
    curated_fusions = [table2cell(data(:, gene_a_col)), ...
                       table2cell(data(:, gene_b_col))];
end

function bed_matrix = load_bed_file(bed_file, chsize, genome_build)
    % Load a BED3-style file with canonical columns chr, start, end.
    % BED format is 0-based, half-open [start, end). Headerless files are also accepted.
    
    validate_file_for_genome_build(bed_file, genome_build);

    data = read_bed_table(bed_file);
    
    bed_matrix = [];
    dropped_rows = 0;
    for i = 1:height(data)
        chr_name = data.chr{i};
        chr_num = parse_chr_name(chr_name);
        start_pos = data.start(i);
        end_pos = data.end(i);
        
        if is_canonical_chr_num(chr_num)
            bed_matrix = [bed_matrix; chr_num, start_pos, end_pos];
        else
            dropped_rows = dropped_rows + 1;
        end
    end

    if dropped_rows > 0
        fprintf('[load_track_files] Ignored %d noncanonical rows in %s\n', dropped_rows, bed_file);
    end
end

function chr_num = parse_chr_name(chr_name)
    % Parse chromosome names to canonical numeric indices:
    % chr1-22/1-22 -> 1-22, chrX/X/23 -> 23, chrY/Y/24 -> 24.
    % Noncanonical chromosomes return -1.
    
    chr_str = lower(string(chr_name));
    chr_str = strtrim(chr_str);
    
    % Remove 'chr' prefix if present
    if startsWith(chr_str, 'chr')
        chr_str = extractAfter(chr_str, 3);
    end
    
    if chr_str == "x" || chr_str == "23"
        chr_num = 23;
    elseif chr_str == "y" || chr_str == "24"
        chr_num = 24;
    else
        chr_num = str2double(chr_str);
        if isnan(chr_num) || chr_num < 1 || chr_num > 24 || floor(chr_num) ~= chr_num
            chr_num = -1;
        end
    end
end

function tf = is_canonical_chr_num(chr_num)
    tf = ~isnan(chr_num) && chr_num >= 1 && chr_num <= 24 && floor(chr_num) == chr_num;
end

function col_idx = find_column_by_alias(tbl, aliases)
    % Find a table column by trying multiple possible names
    
    col_idx = [];
    var_names = lower(string(tbl.Properties.VariableNames));
    
    for i = 1:numel(aliases)
        alias = lower(string(aliases{i}));
        match_idx = find(var_names == alias);
        if ~isempty(match_idx)
            col_idx = match_idx(1);
            return;
        end
    end
end

function col_idx = get_required_column(tbl, aliases, file_path, file_label)
    col_idx = find_column_by_alias(tbl, aliases);
    if isempty(col_idx)
        error('load_track_files:MissingColumns', ...
              'Missing required column in %s (%s). Expected one of: %s', ...
              file_label, file_path, strjoin(aliases, ', '));
    end
end

function data = read_chr_size_table(file_path)
    if file_has_header(file_path, {'chr', 'chrom', 'chromosome', 'size', 'chrom_size', 'length'})
        data = readtable(file_path, 'FileType', 'text', 'Delimiter', '\t', ...
                         'ReadVariableNames', true);
        chr_col = get_required_column(data, {'chr', 'chrom', 'chromosome'}, file_path, 'chrom_sizes_file');
        size_col = get_required_column(data, {'size', 'chrom_size', 'length'}, file_path, 'chrom_sizes_file');
        data = data(:, [chr_col, size_col]);
        data.Properties.VariableNames = {'chr', 'size'};
    else
        data = readtable(file_path, 'FileType', 'text', 'Delimiter', '\t', ...
                         'VariableNames', {'chr', 'size'}, 'ReadVariableNames', false);
    end
end

function data = read_bed_table(file_path)
    if file_has_header(file_path, {'chr', 'chrom', 'chromosome', 'start', 'chromStart', 'end', 'chromEnd'})
        data = readtable(file_path, 'FileType', 'text', 'Delimiter', '\t', ...
                         'ReadVariableNames', true);
        chr_col = get_required_column(data, {'chr', 'chrom', 'chromosome'}, file_path, 'BED file');
        start_col = get_required_column(data, {'start', 'chromStart'}, file_path, 'BED file');
        end_col = get_required_column(data, {'end', 'chromEnd'}, file_path, 'BED file');
        data = data(:, [chr_col, start_col, end_col]);
        data.Properties.VariableNames = {'chr', 'start', 'end'};
    else
        data = readtable(file_path, 'FileType', 'text', 'Delimiter', '\t', ...
                         'VariableNames', {'chr', 'start', 'end'}, ...
                         'ReadVariableNames', false, 'HeaderLines', 0);
    end
end

function tf = file_has_header(file_path, known_header_tokens)
    fid = fopen(file_path, 'r');
    if fid == -1
        error('load_track_files:OpenFailed', 'Could not open file: %s', file_path);
    end

    cleaner = onCleanup(@() fclose(fid));
    first_line = fgetl(fid);
    if ~ischar(first_line)
        tf = false;
        return;
    end

    tokens = split(string(strtrim(first_line)), '\t');
    normalized_tokens = lower(strtrim(tokens));
    normalized_headers = lower(string(known_header_tokens));
    tf = any(ismember(normalized_tokens, normalized_headers));
end

function validate_file_for_genome_build(file_path, genome_build)
    path_lower = lower(file_path);
    build_lower = lower(genome_build);

    if ~contains(path_lower, build_lower)
        fprintf('[load_track_files] WARNING: File path does not contain genome build (%s): %s\n', ...
                genome_build, file_path);
        fprintf('[load_track_files]   Proceeding anyway. Verify file matches expected build.\n');
    end
end
