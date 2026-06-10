%%% Author: Shu Zhang, shu@broadinstitute.org, shu.zhang@gladstone.ucsf.edu
%%% Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu
%%% Date last updated: April 1, 2026

%%% License: GNU AGPL3, Copyright (C) 2023 Dana-Farber Cancer Institute
%%% Dependencies:
%%% See README for guide on how to run this package

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run2DModel
%
% Entry point for running SVsig 2D model.
% Supports both:
% 1) MATLAB struct input for interactive use.
% 2) CLI flags for compiled usage (mcc -m run2DModel.m).
%
% Usage examples:
%   cfg = struct();
%   cfg.work_dir = '/path/to/SVsig';
%   cfg.sv_file = '/path/to/input.csv';
%   cfg.output_file = '/path/to/output.tsv';
%   cfg.model_exist = false;  % set true only when cfg.model_file is provided
%   run2DModel(cfg);
%
%   % CLI:
%   % ./run2DModel -wd /path/to/SVsig -sv /path/to/input.csv -out /path/to/output.tsv
%   % ./run2DModel -help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hits_table = run2DModel(varargin)

    if nargin < 1
        cfg = struct();
    elseif nargin == 1 && isstruct(varargin{1})
        cfg = varargin{1};
    else
        cfg = parseCliArgs(varargin);
        if getOpt(cfg, 'show_help', false)
            printHelp();
            hits_table = table();
            return;
        end
    end

    % Resolve project root relative to this file so there are no hardcoded
    % machine-specific absolute paths.
    this_file_dir = fileparts(mfilename('fullpath'));
    repo_root_dir = fileparts(this_file_dir);

    % I/O settings
    work_dir = getOpt(cfg, 'work_dir', repo_root_dir);
    sv_file = getOpt(cfg, 'sv_file', '');
    output_file = getOpt(cfg, 'output_file', '');

    if ~isfolder(work_dir)
        error('run2DModel:InvalidWorkDir', 'work_dir does not exist: %s', work_dir);
    end
    if isempty(sv_file)
        error('run2DModel:MissingInputFile', 'cfg.sv_file is required and cannot be empty.');
    end
    if ~isfile(sv_file)
        error('run2DModel:MissingInputFile', 'sv_file does not exist: %s', sv_file);
    end
    if isempty(output_file)
        error('run2DModel:MissingOutputFile', 'cfg.output_file is required and cannot be empty.');
    end

    % Model settings
    model_exist = getOpt(cfg, 'model_exist', false);
    model_file = getOpt(cfg, 'model_file', '');  % Required when model_exist = true
    complex_model = getOpt(cfg, 'complex_model', false);
    weights = getOpt(cfg, 'weights', false);
    if complex_model
        weights = true;
    end

    % Filtering and statistical settings
    std_filter = getOpt(cfg, 'std_filter', 10);
    len_filter = getOpt(cfg, 'len_filter', 1e6);
    fdr_threshold = getOpt(cfg, 'fdr_threshold', 0.1);
    tier_std_cutoff = getOpt(cfg, 'tier_std_cutoff', 42833.6);
    low_density_threshold_value = getOpt(cfg, 'low_density_threshold', 5e-5);
    supercluster_distance_threshold = getOpt(cfg, 'supercluster_distance_threshold', 5e4);

    % Track file paths (optional overrides; fall back to data/tracks/<genome_build>/)
    ref_genes_file = getOpt(cfg, 'ref_genes_file', '');
    cancer_genes_file = getOpt(cfg, 'cancer_genes_file', '');
    curated_fusions_file = getOpt(cfg, 'curated_fusions_file', '');
    chrom_sizes_file = getOpt(cfg, 'chrom_sizes_file', '');
    blacklist_file = getOpt(cfg, 'blacklist_file', '');
    l1_elements_file = getOpt(cfg, 'l1_elements_file', '');

    % Binning and genome settings
    num_breakpoints_per_bin = getOpt(cfg, 'num_breakpoints_per_bin', 100);
    bin_length_value = getOpt(cfg, 'bin_length', 5e5);
    n_binshifts = getOpt(cfg, 'n_binshifts', 0);
    max_workers_override = getOpt(cfg, 'max_workers', []);
    chr_list = getOpt(cfg, 'chr_list', 1:23);
    genome_build_value = getOpt(cfg, 'genome_build', 'hg_19');

    if floor(n_binshifts) ~= n_binshifts || n_binshifts < 0
        error('run2DModel:BadBinShiftCount', 'n_binshifts must be a non-negative integer.');
    end
    if ~isempty(max_workers_override)
        if floor(max_workers_override) ~= max_workers_override || max_workers_override < 1
            error('run2DModel:BadMaxWorkers', 'max_workers must be a positive integer when provided.');
        end
    end
    if ~isfinite(supercluster_distance_threshold) || supercluster_distance_threshold < 0
        error('run2DModel:BadSuperclusterDistance', 'supercluster_distance_threshold must be >= 0.');
    end

    % Resolve track file paths: user overrides fall back to data/tracks/<genome_build>/
    user_track_overrides = struct( ...
        'ref_genes_file', ref_genes_file, ...
        'cancer_genes_file', cancer_genes_file, ...
        'curated_fusions_file', curated_fusions_file, ...
        'chrom_sizes_file', chrom_sizes_file, ...
        'blacklist_file', blacklist_file, ...
        'l1_elements_file', l1_elements_file);
    resolved_tracks = resolve_track_paths(work_dir, genome_build_value, user_track_overrides);
    fprintf('[run2DModel] resolved track files:\n');
    fprintf('  ref_genes_file=%s\n', resolved_tracks.ref_genes_file);
    fprintf('  cancer_genes_file=%s\n', resolved_tracks.cancer_genes_file);
    fprintf('  curated_fusions_file=%s\n', resolved_tracks.curated_fusions_file);
    fprintf('  chrom_sizes_file=%s\n', resolved_tracks.chrom_sizes_file);
    fprintf('  blacklist_file=%s\n', resolved_tracks.blacklist_file);
    fprintf('  l1_elements_file=%s\n', resolved_tracks.l1_elements_file);

    % Optionally set a fixed random seed for reproducibility.
    random_seed = getOpt(cfg, 'random_seed', []);
    if ~isempty(random_seed)
        rng(random_seed);
    end

    save_bin_index = getOpt(cfg, 'save_bin_index', false);

    % Print parameters to log file for reproducibility and debugging.
    printRunParameters(work_dir, sv_file, output_file, model_exist, model_file, ...
                       complex_model, weights, std_filter, len_filter, ...
                       fdr_threshold, tier_std_cutoff, ...
                       num_breakpoints_per_bin, bin_length_value, n_binshifts, max_workers_override, chr_list, ...
                       genome_build_value, random_seed, save_bin_index, ...
                       low_density_threshold_value, supercluster_distance_threshold, resolved_tracks);

    global WorkDir
    WorkDir = work_dir;

    global complex
    complex = complex_model;

    global FDR_THRESHOLD
    FDR_THRESHOLD = fdr_threshold;

    global LOW_DENSITY_THRESHOLD
    LOW_DENSITY_THRESHOLD = low_density_threshold_value;

    global CHR
    CHR = chr_list;

    global genome_build
    genome_build = genome_build_value;

    global bin_length
    bin_length = bin_length_value;

    % Set globals for track files so break_invasion_model can access them
    global TRACK_PATHS
    TRACK_PATHS = resolved_tracks;

    addpath(genpath(work_dir));

    % Create output directory when needed.
    output_dir = fileparts(output_file);
    if ~isempty(output_dir) && ~isfolder(output_dir)
        mkdir(output_dir);
    end

    run_context = struct( ...
        'work_dir', work_dir, ...
        'chr_list', chr_list, ...
        'genome_build', genome_build_value, ...
        'track_paths', resolved_tracks, ...
        'low_density_threshold', low_density_threshold_value);

    % Build the shift schedule:
    % run 0 is always the baseline (no shift), and additional runs are
    % fractional shifts of bin_length (k/(n_binshifts+1)).
    shift_values = computeShiftValues(bin_length_value, n_binshifts);
    total_runs = numel(shift_values);

    % Decide parallel capacity: detect available CPUs, then optionally cap
    % with user-provided max_workers.
    cpu_available = detectCpuCount();
    target_workers = chooseWorkerCount(total_runs, cpu_available, max_workers_override);

    fprintf('[run2DModel] n_binshifts=%d, total_runs=%d\n', n_binshifts, total_runs);
    fprintf('[run2DModel] cpu_available=%d, target_workers=%d\n', cpu_available, target_workers);

    hits_by_shift = cell(total_runs, 1);
    bins_by_shift = cell(total_runs, 1);

    % Prefer parallel execution when there are multiple runs and enough
    % worker capacity. Fall back to serial on toolbox/license/runtime issues.
    ran_in_parallel = false;
    if total_runs > 1 && target_workers > 1
        [can_parallel, parallel_reason] = canUseParallelToolbox();
        if can_parallel
            try
                ensureParpoolSize(target_workers);
                parfor run_idx = 1:total_runs
                    shift_bp = shift_values(run_idx);
                    % A zero shift means legacy binning mode; nonzero shifts
                    % trigger the shifted-bin path in SetBins_bylength.
                    shift_arg = normalizeShiftArg(shift_bp);
                    [hits_by_shift{run_idx}, bins_by_shift{run_idx}] = runSVsig( ...
                        sv_file, model_exist, complex_model, weights, len_filter, ...
                        fdr_threshold, bin_length_value, num_breakpoints_per_bin, ...
                        std_filter, model_file, tier_std_cutoff, run_context, shift_arg);
                end
                ran_in_parallel = true;
            catch parallel_err
                fprintf('[run2DModel] Parallel execution failed (%s). Falling back to serial execution.\n', parallel_err.message);
            end
        else
            fprintf('[run2DModel] Parallel execution unavailable: %s\n', parallel_reason);
        end
    end

    if ~ran_in_parallel
        % Serial fallback path (or the only path when one run is requested).
        for run_idx = 1:total_runs
            shift_bp = shift_values(run_idx);
            shift_arg = normalizeShiftArg(shift_bp);
            [hits_by_shift{run_idx}, bins_by_shift{run_idx}] = runSVsig( ...
                sv_file, model_exist, complex_model, weights, len_filter, ...
                fdr_threshold, bin_length_value, num_breakpoints_per_bin, ...
                std_filter, model_file, tier_std_cutoff, run_context, shift_arg);
        end
    end

    if total_runs == 1
        % Preserve exact legacy output schema when no additional shifts are requested.
        hits_table = hits_by_shift{1};
        bins = bins_by_shift{1};
    else
        % Merge shifted runs from already-filtered per-run outputs.
        % Final output is one row per unique SV.
        hits_table = mergeShiftedRunsIntoSuperclusters(hits_by_shift, shift_values, ...
            supercluster_distance_threshold, std_filter, tier_std_cutoff);
        bins = bins_by_shift;
    end

    % Ensure consistent column order and schema for both single-run and multi-run.
    hits_table = matchOutputSchema(hits_table);

    % Write output table
    writetable(hits_table, output_file, 'delimiter', '\t');

    % Optionally save bin indices to sidecar files.
    if save_bin_index
        [out_dir, out_base, out_ext] = fileparts(output_file);
        if total_runs == 1
            % Legacy single-run sidecar filename.
            bin_index_file = fullfile(out_dir, [out_base, out_ext, '.bin_indices.txt']);
            writeBinIndicesFile(bins, bin_index_file);
            fprintf('[run2DModel] bin indices written to: %s\n', bin_index_file);
        else
            % Multi-run sidecars include shift metadata in filenames.
            for run_idx = 1:total_runs
                shift_bp = shift_values(run_idx);
                bin_index_file = fullfile(out_dir, [out_base, out_ext, sprintf('.bin_indices.shift%d_bp%d.txt', run_idx - 1, shift_bp)]);
                writeBinIndicesFile(bins_by_shift{run_idx}, bin_index_file);
                fprintf('[run2DModel] bin indices written to: %s\n', bin_index_file);
            end
        end
    end

end

function cfg = parseCliArgs(args)
    cfg = struct();
    i = 1;
    while i <= numel(args)
        key = string(args{i});

        if key == "-h" || key == "--help" || key == "-help"
            cfg.show_help = true;
            i = i + 1;
            continue;
        end

        if i == numel(args)
            error('run2DModel:BadCliArgs', 'Missing value for flag: %s', key);
        end

        value = args{i + 1};
        switch char(key)
            case {'-wd', '--work_dir'}
                cfg.work_dir = char(value);
            case {'-sv', '--sv_file'}
                cfg.sv_file = char(value);
            case {'-out', '--output_file'}
                cfg.output_file = char(value);
            case {'-model_exist', '--model_exist'}
                cfg.model_exist = parseLogical(value);
            case {'-model', '--model_file'}
                cfg.model_file = char(value);
            case {'-complex', '--complex_model'}
                cfg.complex_model = parseLogical(value);
            case {'-weights', '--weights'}
                cfg.weights = parseLogical(value);
            case {'-std', '--std_filter'}
                cfg.std_filter = str2double(char(value));
            case {'-tier', '--tier_std_cutoff'}
                cfg.tier_std_cutoff = str2double(char(value));
            case {'-ldt', '--low_density_threshold'}
                cfg.low_density_threshold = str2double(char(value));
            case {'-scdist', '--supercluster_distance_threshold'}
                cfg.supercluster_distance_threshold = str2double(char(value));
            case {'-len', '--len_filter'}
                cfg.len_filter = str2double(char(value));
            case {'-fdr', '--fdr_threshold'}
                cfg.fdr_threshold = str2double(char(value));
            case {'-nbpb', '--num_breakpoints_per_bin'}
                cfg.num_breakpoints_per_bin = str2double(char(value));
            case {'-bin', '--bin_length'}
                cfg.bin_length = str2double(char(value));
            case {'-nshift', '--n_binshifts'}
                cfg.n_binshifts = str2double(char(value));
            case {'-maxw', '--max_workers'}
                cfg.max_workers = str2double(char(value));
            case {'-chr', '--chr_list'}
                cfg.chr_list = parseChrList(char(value));
            case {'-genome', '--genome_build'}
                cfg.genome_build = char(value);
            case {'-seed', '--random_seed'}
                cfg.random_seed = str2double(char(value));
            case {'-bix', '--save_bin_index'}
                cfg.save_bin_index = parseLogical(value);
            case {'-ref_genes', '--ref_genes_file'}
                cfg.ref_genes_file = char(value);
            case {'-cancer_genes', '--cancer_genes_file'}
                cfg.cancer_genes_file = char(value);
            case {'-curated_fusions', '--curated_fusions_file'}
                cfg.curated_fusions_file = char(value);
            case {'-chrom_sizes', '--chrom_sizes_file'}
                cfg.chrom_sizes_file = char(value);
            case {'-blacklist', '--blacklist_file'}
                cfg.blacklist_file = char(value);
            case {'-l1', '--l1_elements_file'}
                cfg.l1_elements_file = char(value);
            otherwise
                error('run2DModel:BadCliArgs', 'Unknown flag: %s', key);
        end

        i = i + 2;
    end
end

function tf = parseLogical(value)
    if islogical(value)
        tf = value;
        return;
    end

    text = lower(string(value));
    if text == "1" || text == "true" || text == "yes"
        tf = true;
    elseif text == "0" || text == "false" || text == "no"
        tf = false;
    else
        error('run2DModel:BadCliArgs', 'Invalid logical value: %s', string(value));
    end
end

function chr_list = parseChrList(value)
    % Accept numeric chromosome specifications only:
    %   1:23 (range form)
    %   1:24 (range form, optional)
    %   1,2,3,23 (comma-separated form)
    % All values must be integers in [1, 24].
    value = strrep(value, ' ', '');
    
    if contains(value, ':')
        parts = split(value, ':');
        if numel(parts) ~= 2
            error('run2DModel:BadCliArgs', 'Invalid chr_list value: %s', value);
        end
        start_v = str2double(parts{1});
        end_v = str2double(parts{2});
        if isnan(start_v) || isnan(end_v) || start_v < 1 || end_v > 24 || ...
           floor(start_v) ~= start_v || floor(end_v) ~= end_v || start_v > end_v
            error('run2DModel:BadCliArgs', 'Invalid chr_list value: %s', value);
        end
        chr_list = start_v:end_v;
    else
        parts = split(value, ',');
        chr_list = zeros(1, numel(parts));
        for k = 1:numel(parts)
            chr_num = str2double(parts{k});
            if isnan(chr_num) || chr_num < 1 || chr_num > 24 || floor(chr_num) ~= chr_num
                error('run2DModel:BadCliArgs', 'Invalid chr_list value: %s', value);
            end
            chr_list(k) = chr_num;
        end
    end
end

function printHelp()
    fprintf('\nUsage: run2DModel [OPTIONS]\n\n');
    fprintf('Required:\n');
    fprintf('  -sv,  --sv_file                    [path]   Input SV CSV file.\n');
    fprintf('  -out, --output_file                [path]   Output TSV file to write.\n\n');

    fprintf('Model Parameters:\n');
    fprintf('  -model_exist, --model_exist        [bool]   If true, load a pre-computed model_file.\n');
    fprintf('                                              If false, compute background model from input.\n');
    fprintf('                                              default=false\n');
    fprintf('  -model, --model_file               [path]   Pre-computed background model (.mat).\n');
    fprintf('                                              Required only when model_exist=true.\n\n');
    fprintf('  -complex, --complex_model          [bool]   Enable complex model mode.\n');
    fprintf('                                              default=false\n');
    fprintf('  -weights, --weights                [bool]   Use weights for individual rearrangements.\n');
    fprintf('                                              default=false (forced true when complex_model=true)\n\n');

    fprintf('Binning and Clustering:\n');
    fprintf('  -bin, --bin_length                 [num]    Length of bins for adjacency and binning.\n');
    fprintf('                                              default=5e5\n');
    fprintf('  -nshift, --n_binshifts             [int]    Number of shifted bin runs to add.\n');
    fprintf('                                              0 => no shift (single run),\n');
    fprintf('                                              1 => add 1/2-bin shift (2 total runs),\n');
    fprintf('                                              2 => add 1/3 and 2/3 shifts (3 total runs), etc.\n');
    fprintf('                                              default=0\n');
    fprintf('  -nbpb, --num_breakpoints_per_bin   [num]    Minimum breakpoints per bin for adjacency and binning.\n');
    fprintf('                                              default=100\n\n');


    fprintf('Filtering and Significance:\n');
    fprintf('  -std, --std_filter                 [num]    Keep hits only if stddev on both breakpoints exceeds this.\n');
    fprintf('                                              default=10\n');
    fprintf('  -len, --len_filter                 [num]    Minimum span used in downstream hit filtering.\n');
    fprintf('                                              default=1e6\n');
    fprintf('  -fdr, --fdr_threshold              [num]    BH-FDR q-value threshold for significant tiles.\n');
    fprintf('                                              default=0.1\n');
    fprintf('  -tier, --tier_std_cutoff           [num]    Cutoff used to assign output tier1/2/3 from stddev_i and stddev_j.\n');
    fprintf('                                              default=42833.6\n\n');
    fprintf('  -ldt, --low_density_threshold      [num]    Threshold used to remove low-density bins.\n');
    fprintf('                                              default=5e-5\n\n');
    fprintf('  -scdist, --supercluster_distance_threshold [num]  Max breakpoint distance (bp) used to\n');
    fprintf('                                              merge clusters across shifted runs into superclusters.\n');
    fprintf('                                              default=50000\n\n');

    fprintf('Track Files (optional; fall back to data/tracks/<genome_build>/):\n');
    fprintf('  -ref_genes, --ref_genes_file       [path]   Reference genes file (default: ref_genes.txt.gz)\n');
    fprintf('  -cancer_genes, --cancer_genes_file [path]   Cancer genes file (default: cancer_genes.tsv)\n');
    fprintf('  -curated_fusions, --curated_fusions_file [path]   Curated fusions file (default: curated_fusions.tsv)\n');
    fprintf('  -chrom_sizes, --chrom_sizes_file   [path]   Chromosome sizes file (default: chrom.sizes)\n');
    fprintf('  -blacklist, --blacklist_file       [path]   Blacklist regions file (default: blacklist.bed)\n');
    fprintf('  -l1, --l1_elements_file            [path]   L1 elements file (default: l1_elements.bed)\n\n');

    fprintf('Genome and Runtime:\n');
    fprintf('  -chr, --chr_list                   [list]   Chromosomes included in model building and testing.\n');
    fprintf('                                              Accepts numeric form: 1:23, 1:24, or 1,2,3,...,23\n');
    fprintf('                                              default=1:23\n');
    fprintf('  -genome, --genome_build            [text]   Genome build used by downstream annotation.\n');
    fprintf('                                              default=hg_19\n');
    fprintf('  -wd, --work_dir                    [path]   Root directory added to MATLAB path.\n');
    fprintf('                                              default=<repo root inferred from run2DModel.m>\n');
    fprintf('  -seed, --random_seed               [int]    Optional RNG seed for reproducibility.\n');
    fprintf('                                              default=unset\n');
    fprintf('  -maxw, --max_workers               [int]    Maximum parallel workers for shifted runs.\n');
    fprintf('                                              default=auto (infer from available CPUs)\n');
    fprintf('  -bix, --save_bin_index             [bool]   Write a sidecar bin index file per shift alongside the output.\n');
    fprintf('                                              Contains chr, start, and end for each unique bin in the\n');
    fprintf('                                              final filtered output table for that shift.\n');
    fprintf('                                              Filename: <output>.bin_indices.shift<N>bp<ext>.\n');
    fprintf('                                              default=false\n');
    fprintf('  -h, --help                                  Show this help and exit.\n\n');
end

function value = getOpt(cfg, field_name, default_value)
    if isfield(cfg, field_name)
        value = cfg.(field_name);
    else
        value = default_value;
    end
end

function printRunParameters(work_dir, sv_file, output_file, model_exist, model_file, ...
                            complex_model, weights, std_filter, len_filter, ...
                            fdr_threshold, tier_std_cutoff, ...
                            num_breakpoints_per_bin, bin_length_value, n_binshifts, max_workers_override, chr_list, ...
                            genome_build_value, random_seed, save_bin_index, ...
                            low_density_threshold_value, supercluster_distance_threshold, resolved_tracks)
    fprintf('\n[run2DModel] effective parameters\n');
    fprintf('--work_dir=%s\n', work_dir);
    fprintf('--sv_file=%s\n', sv_file);
    fprintf('--output_file=%s\n', output_file);
    fprintf('--model_exist=%s\n', boolToText(model_exist));
    fprintf('--model_file=%s\n', emptyToNA(model_file));
    fprintf('--complex_model=%s\n', boolToText(complex_model));
    fprintf('--weights=%s\n', boolToText(weights));
    fprintf('--std_filter=%g\n', std_filter);
    fprintf('--len_filter=%g\n', len_filter);
    fprintf('--fdr_threshold=%g\n', fdr_threshold);
    fprintf('--tier_std_cutoff=%g\n', tier_std_cutoff);
    fprintf('--low_density_threshold=%g\n', low_density_threshold_value);
    fprintf('--supercluster_distance_threshold=%g\n', supercluster_distance_threshold);
    fprintf('--num_breakpoints_per_bin=%g\n', num_breakpoints_per_bin);
    fprintf('--bin_length=%g\n', bin_length_value);
    fprintf('--n_binshifts=%g\n', n_binshifts);
    fprintf('--max_workers=%s\n', maxWorkersToText(max_workers_override));
    fprintf('--chr_list=%s\n', chrListToText(chr_list));
    fprintf('--genome_build=%s\n', genome_build_value);
    fprintf('--random_seed=%s\n', seedToText(random_seed));
    fprintf('--save_bin_index=%s\n', boolToText(save_bin_index));
    fprintf('--resolved track files:\n');
    fprintf('  ref_genes_file=%s\n', emptyToNA(resolved_tracks.ref_genes_file));
    fprintf('  cancer_genes_file=%s\n', emptyToNA(resolved_tracks.cancer_genes_file));
    fprintf('  curated_fusions_file=%s\n', emptyToNA(resolved_tracks.curated_fusions_file));
    fprintf('  chrom_sizes_file=%s\n', emptyToNA(resolved_tracks.chrom_sizes_file));
    fprintf('  blacklist_file=%s\n', emptyToNA(resolved_tracks.blacklist_file));
    fprintf('  l1_elements_file=%s\n\n', emptyToNA(resolved_tracks.l1_elements_file));
end

function text = boolToText(tf)
    if tf
        text = 'true';
    else
        text = 'false';
    end
end

function text = emptyToNA(value)
    if isempty(value)
        text = 'NA';
    else
        text = value;
    end
end

function text = seedToText(random_seed)
    if isempty(random_seed)
        text = 'NA';
    else
        text = num2str(random_seed);
    end
end

function text = maxWorkersToText(max_workers_override)
    if isempty(max_workers_override)
        text = 'auto';
    else
        text = num2str(max_workers_override);
    end
end

function text = chrListToText(chr_list)
    if isempty(chr_list)
        text = 'NA';
        return;
    end

    if isvector(chr_list) && numel(chr_list) > 1 && all(diff(chr_list) == 1)
        text = sprintf('%g:%g', chr_list(1), chr_list(end));
    else
        text = sprintf('%g,', chr_list);
        text = text(1:end-1);
    end
end

function shift_values = computeShiftValues(bin_length_value, n_binshifts)
    % Returns [0, floor(1/(n+1)*bin), ..., floor(n/(n+1)*bin)].
    % The first value (0) is the baseline no-shift run.
    shift_values = zeros(1, n_binshifts + 1);
    for idx = 1:(n_binshifts + 1)
        shift_values(idx) = floor(bin_length_value * (idx - 1) / (n_binshifts + 1));
    end
end

function shift_arg = normalizeShiftArg(shift_bp)
    % Keep baseline run identical to legacy behavior by passing [] for 0.
    if shift_bp <= 0
        shift_arg = [];
    else
        shift_arg = shift_bp;
    end
end

function cpu_count = detectCpuCount()
    % Prefer MATLAB core detection; fall back to JVM processor count.
    cpu_count = 1;
    try
        if exist('feature', 'builtin') == 5
            cpu_count = feature('numcores');
        end
    catch
    end

    if isempty(cpu_count) || ~isfinite(cpu_count) || cpu_count < 1
        try
            cpu_count = java.lang.Runtime.getRuntime().availableProcessors();
        catch
            cpu_count = 1;
        end
    end

    cpu_count = max(1, floor(double(cpu_count)));
end

function worker_count = chooseWorkerCount(total_runs, cpu_available, max_workers_override)
    % Never request more workers than runs, detected CPUs, or user cap.
    worker_count = min(total_runs, cpu_available);
    if ~isempty(max_workers_override)
        worker_count = min(worker_count, max_workers_override);
    end
    worker_count = max(1, floor(worker_count));
end

function [can_parallel, reason] = canUseParallelToolbox()
    % Validate both runtime support and license before using parfor/parpool.
    can_parallel = true;
    reason = '';

    if exist('parfor', 'builtin') ~= 5 && exist('parfor', 'file') ~= 2
        can_parallel = false;
        reason = 'parfor is unavailable in this MATLAB environment';
        return;
    end

    if license('test', 'Distrib_Computing_Toolbox') == 0
        can_parallel = false;
        reason = 'Parallel Computing Toolbox license unavailable';
    end
end

function ensureParpoolSize(target_workers)
    % Reuse existing pool when size matches; otherwise recreate at target size.
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(target_workers);
        return;
    end

    if pool.NumWorkers ~= target_workers
        delete(pool);
        parpool(target_workers);
    end
end

function merged_hits = mergeShiftedRunsIntoSuperclusters(hits_by_shift, shift_values, ...
    distance_threshold, std_filter, tier_std_cutoff)
    % Stack all shifted runs and tag each row with run-specific metadata.
    all_hits = table();
    for run_idx = 1:numel(hits_by_shift)
        run_table = hits_by_shift{run_idx};
        if isempty(run_table)
            continue;
        end

        shift_run_value = run_idx - 1;
        run_table.shift_run = repmat(shift_run_value, height(run_table), 1);
        run_table.bin_shift_bp = repmat(shift_values(run_idx), height(run_table), 1);
        run_table.cluster_label = cellstr("s" + string(shift_run_value) + "_c" + string(run_table.cluster_num));
        all_hits = [all_hits; run_table]; %#ok<AGROW>
    end

    if isempty(all_hits)
        merged_hits = all_hits;
        return;
    end

    % Build a cluster-level graph from row proximity across shifted outputs.
    [cluster_labels, ~, row_cluster_idx] = unique(all_hits.cluster_label, 'stable');
    edge_pairs = buildSuperclusterEdges(all_hits, row_cluster_idx, distance_threshold);

    num_clusters = numel(cluster_labels);
    if isempty(edge_pairs)
        cluster_component = (1:num_clusters)';
    else
        g = graph(edge_pairs(:,1), edge_pairs(:,2), [], num_clusters);
        cluster_component = conncomp(g)';
    end

    row_component = cluster_component(row_cluster_idx);
    num_components = max(row_component);

    shifts_found_by_component = zeros(num_components, 1);
    for comp_idx = 1:num_components
        rows_in_component = (row_component == comp_idx);
        shifts_found_by_component(comp_idx) = numel(unique(all_hits.shift_run(rows_in_component)));
    end

    all_hits.run_cluster_num = all_hits.cluster_num;
    all_hits.supercluster_id = row_component;
    all_hits.cluster_num = row_component;
    all_hits.shifts_found_count = shifts_found_by_component(row_component);

    % Recompute component-level summary metrics and tier assignments.
    super_nsamples = zeros(num_components, 1);
    super_std_i = zeros(num_components, 1);
    super_std_j = zeros(num_components, 1);
    super_pval = zeros(num_components, 1);
    super_prob = zeros(num_components, 1);
    super_qval = zeros(num_components, 1);
    super_tier_rank = 3 * ones(num_components, 1);
    super_keep = false(num_components, 1);

    for comp_idx = 1:num_components
        rows_in_component = (row_component == comp_idx);
        sub = all_hits(rows_in_component, :);

        super_nsamples(comp_idx) = numel(unique(string(sub.sid)));
        super_std_i(comp_idx) = std(sub.pos_i);
        super_std_j(comp_idx) = std(sub.pos_j);
        super_pval(comp_idx) = min(sub.pval);
        super_qval(comp_idx) = min(sub.tile_qval);
        super_prob(comp_idx) = mean(sub.prob, 'omitnan');

        if super_std_i(comp_idx) <= tier_std_cutoff && super_std_j(comp_idx) <= tier_std_cutoff
            super_tier_rank(comp_idx) = 1;
        elseif xor(super_std_i(comp_idx) <= tier_std_cutoff, super_std_j(comp_idx) <= tier_std_cutoff)
            super_tier_rank(comp_idx) = 2;
        else
            super_tier_rank(comp_idx) = 3;
        end

        super_keep(comp_idx) = (super_nsamples(comp_idx) > 1) && ...
                              (super_std_i(comp_idx) > std_filter) && ...
                              (super_std_j(comp_idx) > std_filter);
    end

    all_hits.num_hits = super_nsamples(row_component);
    all_hits.stddev_i = super_std_i(row_component);
    all_hits.stddev_j = super_std_j(row_component);
    all_hits.pval = super_pval(row_component);
    all_hits.tile_qval = super_qval(row_component);
    all_hits.prob = super_prob(row_component);
    all_hits.supercluster_nsamples = super_nsamples(row_component);
    all_hits.supercluster_pval = super_pval(row_component);
    all_hits.supercluster_prob = super_prob(row_component);
    all_hits.supercluster_tier = cellstr(mapTierRankToText(super_tier_rank(row_component)));

    % Track where each unique SV appeared across runs/clusters so the
    % per-run output can be reconstructed downstream.
    sv_key_cols = {'sid','chr_i','pos_i','strand_i','chr_j','pos_j','strand_j'};
    [unique_sv_keys, ~, sv_idx] = unique(all_hits(:, sv_key_cols), 'rows', 'stable');
    sv_cluster_ids = strings(height(unique_sv_keys), 1);
    sv_cluster_tiers = strings(height(unique_sv_keys), 1);
    sv_tile_qval = strings(height(unique_sv_keys), 1);
    sv_pval = strings(height(unique_sv_keys), 1);
    sv_prob = strings(height(unique_sv_keys), 1);
    supercluster_cluster_tuples = strings(num_components, 1);

    % Compute supercluster-level cluster tuples and per-SV tier mapping.
    for comp_idx = 1:num_components
        rows_in_component = (row_component == comp_idx);
        tuple_rows = [all_hits.shift_run(rows_in_component), all_hits.run_cluster_num(rows_in_component)];
        tuple_rows = unique(tuple_rows, 'rows', 'stable');
        tuple_text = arrayfun(@(r) sprintf('(%d, %d)', tuple_rows(r,1), tuple_rows(r,2)), ...
                              1:size(tuple_rows,1), 'UniformOutput', false);
        supercluster_cluster_tuples(comp_idx) = strjoin(tuple_text, ', ');
    end

    for sv_i = 1:height(unique_sv_keys)
        sv_rows = (sv_idx == sv_i);
        tuple_rows = [all_hits.shift_run(sv_rows), all_hits.run_cluster_num(sv_rows)];
        tuple_rows = unique(tuple_rows, 'rows', 'stable');
        tuple_text = arrayfun(@(r) sprintf('(%d, %d)', tuple_rows(r,1), tuple_rows(r,2)), ...
                              1:size(tuple_rows,1), 'UniformOutput', false);

        % Extract tier, tile_qval, pval, prob for each tuple in order, matching the tuple creation.
        sv_tiers_ordered = {};
        sv_tile_qval_ordered = {};
        sv_pval_ordered = {};
        sv_prob_ordered = {};
        for t = 1:size(tuple_rows, 1)
            shift_run_val = tuple_rows(t, 1);
            cluster_num_val = tuple_rows(t, 2);
            idx_match = find(all_hits.shift_run(sv_rows) == shift_run_val & ...
                            all_hits.run_cluster_num(sv_rows) == cluster_num_val, 1, 'first');
            sv_rows_indices = find(sv_rows);
            sv_tiers_ordered{t} = char(all_hits.supercluster_tier{sv_rows_indices(idx_match)});
            sv_tile_qval_ordered{t} = sprintf('%.4g', all_hits.tile_qval(sv_rows_indices(idx_match)));
            sv_pval_ordered{t} = sprintf('%.4g', all_hits.pval(sv_rows_indices(idx_match)));
            sv_prob_ordered{t} = sprintf('%.4g', all_hits.prob(sv_rows_indices(idx_match)));
        end

        sv_cluster_ids(sv_i) = strjoin(tuple_text, ', ');
        sv_cluster_tiers(sv_i) = strjoin(sv_tiers_ordered, ', ');
        sv_tile_qval(sv_i) = strjoin(sv_tile_qval_ordered, ', ');
        sv_pval(sv_i) = strjoin(sv_pval_ordered, ', ');
        sv_prob(sv_i) = strjoin(sv_prob_ordered, ', ');
    end

    % Keep one row per unique SV event identity after supercluster mapping.
    [~, unique_row_idx] = unique(all_hits(:, sv_key_cols), 'rows', 'stable');
    merged_hits = all_hits(unique_row_idx, :);

    [~, merged_sv_idx] = ismember(merged_hits(:, sv_key_cols), unique_sv_keys, 'rows');
    merged_hits.sv_cluster_ids = cellstr(sv_cluster_ids(merged_sv_idx));
    merged_hits.sv_cluster_tiers = cellstr(sv_cluster_tiers(merged_sv_idx));
    merged_hits.sv_tile_qval = cellstr(sv_tile_qval(merged_sv_idx));
    merged_hits.sv_pval = cellstr(sv_pval(merged_sv_idx));
    merged_hits.sv_prob = cellstr(sv_prob(merged_sv_idx));
    merged_hits.supercluster_cluster_ids = cellstr(supercluster_cluster_tuples(merged_hits.cluster_num));

    keep_rows = super_keep(merged_hits.cluster_num);
    merged_hits = merged_hits(keep_rows, :);

    if ~isempty(merged_hits)
        merged_hits = sortrows(merged_hits, {'chr_i', 'pos_i', 'chr_j', 'pos_j'});
    end

    % Drop transient run-level columns and old tier from the final merged output.
    transient_cols = {'shift_run', 'bin_shift_bp', 'cluster_label', 'run_cluster_num', 'tier', 'sv_run_cluster_ids'};
    keep_cols = ~ismember(merged_hits.Properties.VariableNames, transient_cols);
    merged_hits = merged_hits(:, keep_cols);
end

function tier_text = mapTierRankToText(tier_rank)
    tier_text = repmat("tier 3", numel(tier_rank), 1);
    tier_text(tier_rank <= 1) = "tier 1";
    tier_text(tier_rank == 2) = "tier 2";
end

function edge_pairs = buildSuperclusterEdges(all_hits, row_cluster_idx, distance_threshold)
    % Connect clusters if any pair of rows from matching chromosome-pairs
    % are within distance_threshold on both breakpoints.
    edge_pairs = zeros(0, 2);
    if isempty(all_hits)
        return;
    end

    chr_pair = [all_hits.chr_i, all_hits.chr_j];
    [unique_pairs, ~, pair_idx] = unique(chr_pair, 'rows', 'stable'); %#ok<ASGLU>

    for p = 1:size(unique_pairs, 1)
        rows = find(pair_idx == p);
        if numel(rows) < 2
            continue;
        end

        pos_i = all_hits.pos_i(rows);
        pos_j = all_hits.pos_j(rows);
        [pos_i_sorted, order] = sort(pos_i);
        rows_sorted = rows(order);
        pos_j_sorted = pos_j(order);

        n = numel(rows_sorted);
        for a = 1:(n - 1)
            b = a + 1;
            while b <= n && (pos_i_sorted(b) - pos_i_sorted(a)) <= distance_threshold
                if abs(pos_j_sorted(b) - pos_j_sorted(a)) <= distance_threshold
                    c1 = row_cluster_idx(rows_sorted(a));
                    c2 = row_cluster_idx(rows_sorted(b));
                    if c1 ~= c2
                        edge_pairs(end + 1, :) = [min(c1, c2), max(c1, c2)]; %#ok<AGROW>
                    end
                end
                b = b + 1;
            end
        end
    end

    if ~isempty(edge_pairs)
        edge_pairs = unique(edge_pairs, 'rows');
    end
end

function output_table = matchOutputSchema(input_table)
    % Normalize schema so single-run and multi-run outputs have identical structure.
    output_table = input_table;

    % Ensure critical columns exist in expected order.
    expected_cols = {'supercluster_id', 'sv_cluster_ids', 'sv_cluster_tiers', ...
                     'sv_tile_qval', 'sv_pval', 'sv_prob', ...
                     'supercluster_cluster_ids', 'supercluster_tier', 'supercluster_pval', ...
                     'supercluster_prob', 'supercluster_nsamples', 'shifts_found_count'};
    for col = expected_cols
        col_name = col{1};
        if ~any(strcmp(output_table.Properties.VariableNames, col_name))
            output_table.(col_name) = repmat("", height(output_table), 1);
        end
    end

    % Remove redundant or deprecated columns.
    deprecated_cols = {'tier', 'sv_run_cluster_ids', 'sv_shift_cluster_tuples', 'cluster_num', 'num_hits'};
    for col = deprecated_cols
        col_name = col{1};
        if any(strcmp(output_table.Properties.VariableNames, col_name))
            output_table.(col_name) = [];
        end
    end
end

function writeBinIndicesFile(bins, bin_index_file)
    if istable(bins)
        writetable(bins, bin_index_file, 'Delimiter', '\t');
        return;
    end

    if isnumeric(bins) || islogical(bins)
        writematrix(bins, bin_index_file, 'Delimiter', '\t');
        return;
    end

    error('run2DModel:BadBinIndexType', ...
          'Unsupported bins type for sidecar output: %s', class(bins));
end
