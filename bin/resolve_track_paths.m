function track_paths = resolve_track_paths(work_dir, genome_build, user_overrides)
    % Resolve genomic track file paths with user override support.
    % Falls back to standard directory structure if user paths are empty or invalid.
    %
    % Args:
    %   work_dir: SVsig project root
    %   genome_build: 'hg_19' or 'hg_38'
    %   user_overrides: struct with optional fields:
    %       ref_genes_file, cancer_genes_file, curated_fusions_file,
    %       chrom_sizes_file, blacklist_file, l1_elements_file
    %
    % Returns:
    %   track_paths: struct with resolved full paths for each track
    
    % Standard directory for default files
    default_dir = fullfile(work_dir, 'data', 'tracks', genome_build);
    
    % Default filenames (edit these if you change the standard naming)
    defaults = struct(...
        'ref_genes_file', 'hg19.ncbiRefSeq.colnames.txt', ...
        'cancer_genes_file', 'Cosmic_CancerGeneCensus_v104_GRCh37.tsv', ...
        'curated_fusions_file', 'CuratedFusionGene.colnames.txt', ...
        'chrom_sizes_file', 'hg19.chrom.sizes', ...
        'blacklist_file', 'mask_track.bed', ...
        'l1_elements_file', 'l1_track.bed');
    
    % Field names to process
    field_names = fieldnames(defaults);
    track_paths = struct();
    
    for i = 1:numel(field_names)
        field = field_names{i};
        
        % Check if user provided an override
        if isfield(user_overrides, field) && ~isempty(user_overrides.(field))
            user_path = user_overrides.(field);
            if isfile(user_path)
                track_paths.(field) = user_path;
                fprintf('[resolve_track_paths] Using user override for %s: %s\n', field, user_path);
            else
                fprintf('[resolve_track_paths] WARNING: User-supplied %s does not exist: %s\n', field, user_path);
                fprintf('[resolve_track_paths]   Falling back to default: %s\n', fullfile(default_dir, defaults.(field)));
                track_paths.(field) = fullfile(default_dir, defaults.(field));
            end
        else
            % Use default
            track_paths.(field) = fullfile(default_dir, defaults.(field));
        end
    end
    
end
