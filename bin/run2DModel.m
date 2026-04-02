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
    bks_cluster = getOpt(cfg, 'bks_cluster', 1);
    fdr_threshold = getOpt(cfg, 'fdr_threshold', 0.1);
    tier_std_cutoff = getOpt(cfg, 'tier_std_cutoff', 42833.6);

    % Binning and genome settings
    num_breakpoints_per_bin = getOpt(cfg, 'num_breakpoints_per_bin', 100);
    bin_length_value = getOpt(cfg, 'bin_length', 5e5);
    chr_list = getOpt(cfg, 'chr_list', 1:23);
    genome_build_value = getOpt(cfg, 'genome_build', 'hg_19');

    % Optionally set a fixed random seed for reproducibility.
    random_seed = getOpt(cfg, 'random_seed', []);
    if ~isempty(random_seed)
        rng(random_seed);
    end
    
    % Print parameters to log file for reproducibility and debugging.
    printRunParameters(work_dir, sv_file, output_file, model_exist, model_file, ...
                       complex_model, weights, std_filter, len_filter, ...
                       bks_cluster, fdr_threshold, tier_std_cutoff, ...
                       num_breakpoints_per_bin, bin_length_value, chr_list, ...
                       genome_build_value, random_seed);

    global WorkDir
    WorkDir = work_dir;

    global complex
    complex = complex_model;

    global FDR_THRESHOLD
    FDR_THRESHOLD = fdr_threshold;

    global CHR
    CHR = chr_list;

    global genome_build
    genome_build = genome_build_value;

    global bin_length
    bin_length = bin_length_value;

    addpath(genpath(work_dir));

    % Create output directory when needed.
    output_dir = fileparts(output_file);
    if ~isempty(output_dir) && ~isfolder(output_dir)
        mkdir(output_dir);
    end

    % Run model
    hits_table = runSVsig(sv_file, model_exist, complex_model, weights, len_filter, bks_cluster, ...
                          fdr_threshold, bin_length_value, num_breakpoints_per_bin, ...
                          std_filter, model_file, tier_std_cutoff);

    % Write output table
    writetable(hits_table, output_file, 'delimiter', '\t');

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
            case {'-len', '--len_filter'}
                cfg.len_filter = str2double(char(value));
            case {'-bks', '--bks_cluster'}
                cfg.bks_cluster = str2double(char(value));
            case {'-fdr', '--fdr_threshold'}
                cfg.fdr_threshold = str2double(char(value));
            case {'-nbpb', '--num_breakpoints_per_bin'}
                cfg.num_breakpoints_per_bin = str2double(char(value));
            case {'-bin', '--bin_length'}
                cfg.bin_length = str2double(char(value));
            case {'-chr', '--chr_list'}
                cfg.chr_list = parseChrList(char(value));
            case {'-genome', '--genome_build'}
                cfg.genome_build = char(value);
            case {'-seed', '--random_seed'}
                cfg.random_seed = str2double(char(value));
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
    % Accept either "1:23" or comma-separated values like "1,2,3,23".
    if contains(value, ':')
        parts = split(value, ':');
        if numel(parts) ~= 2
            error('run2DModel:BadCliArgs', 'Invalid chr_list value: %s', value);
        end
        start_v = str2double(parts{1});
        end_v = str2double(parts{2});
        chr_list = start_v:end_v;
    else
        parts = split(value, ',');
        chr_list = zeros(1, numel(parts));
        for k = 1:numel(parts)
            chr_list(k) = str2double(parts{k});
        end
    end

    if any(isnan(chr_list))
        error('run2DModel:BadCliArgs', 'Invalid chr_list value: %s', value);
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
    fprintf('  -bks, --bks_cluster                [0|1]    Select breakpoint significance method variant (PVal(fragile sites not accounted for) or PvalMH(fragile sites account for)).\n');
    fprintf('                                              default=1\n');
    fprintf('  -bin, --bin_length                 [num]    Length of bins for adjacency and binning.\n');
    fprintf('                                              default=5e5\n');
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

    fprintf('Genome and Runtime:\n');
    fprintf('  -chr, --chr_list                   [list]   Chromosomes included in model building and testing.\n');
    fprintf('                                              default=1:23\n');
    fprintf('  -genome, --genome_build            [text]   Genome build used by downstream annotation.\n');
    fprintf('                                              default=hg_19\n');
    fprintf('  -wd, --work_dir                    [path]   Root directory added to MATLAB path.\n');
    fprintf('                                              default=<repo root inferred from run2DModel.m>\n');
    fprintf('  -seed, --random_seed               [int]    Optional RNG seed for reproducibility.\n');
    fprintf('                                              default=unset\n');
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
                            bks_cluster, fdr_threshold, tier_std_cutoff, ...
                            num_breakpoints_per_bin, bin_length_value, chr_list, ...
                            genome_build_value, random_seed)
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
    fprintf('--bks_cluster=%g\n', bks_cluster);
    fprintf('--fdr_threshold=%g\n', fdr_threshold);
    fprintf('--tier_std_cutoff=%g\n', tier_std_cutoff);
    fprintf('--num_breakpoints_per_bin=%g\n', num_breakpoints_per_bin);
    fprintf('--bin_length=%g\n', bin_length_value);
    fprintf('--chr_list=%s\n', chrListToText(chr_list));
    fprintf('--genome_build=%s\n', genome_build_value);
    fprintf('--random_seed=%s\n\n', seedToText(random_seed));
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
