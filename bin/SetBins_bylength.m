function [bins, numbins] = SetBins_bylength(events,avg_bin_length,chsize,CHR,bin_shift)

disp('setting bins boundaries...');

chr = [events(:,1);events(:,4)];
pos = [events(:,2);events(:,5)];

bc=1;
bins = []; 

% Parse bin_shift parameter
if nargin >= 5 && ~isempty(bin_shift)
    shift = floor(double(bin_shift));
    if ~isfinite(shift) || shift < 0
        error('SetBins_bylength:BadBinShift', 'bin_shift must be a non-negative integer.');
    end
else
    shift = 0;
end

% Loop through selected chromosomes
for c1 = CHR
    chri = find(chr==c1);
    sorted_pos = sort(pos(chri));
    chr_size = floor(chsize(c1));

    % Pre-calculate the UNIFORM geometry based on the full chromosome
    % This ensures the grid matches the nshift=0 baseline exactly.
    num_chr = floor(chr_size / avg_bin_length);
    if num_chr > 0
        avg_length_chr = floor(avg_bin_length + ((chr_size - num_chr*avg_bin_length)/num_chr));
    else
        num_chr = 0;
        avg_length_chr = chr_size;
    end

    % Apply shift offset padding bin if requested
    if shift > 0
        bins(bc,1) = c1;
        bins(bc,2) = 1;
        bins(bc,3) = min(shift, chr_size);
        
        breaks_per_bin = find(sorted_pos(:,1) <= (bins(bc,3) + 1) & sorted_pos(:,1) >= 1);
        bins(bc,4) = length(breaks_per_bin);
        bc = bc + 1;
        
        start1 = shift + 1;
    else
        start1 = 1;
    end

    % Generate bins using the pre-calculated geometry
    if num_chr > 0
        for c2 = 1:num_chr
            bins(bc,1) = c1;
            bins(bc,2) = start1;
            % Cap at chr_size to maintain structural integrity
            bins(bc,3) = min(start1 + avg_length_chr - 1, chr_size);
            
            breaks_per_bin = find(sorted_pos(:,1) <= (bins(bc,3) + 1) & sorted_pos(:,1) >= start1);
            bins(bc,4) = length(breaks_per_bin);
            bc = bc + 1;
            
            start1 = start1 + avg_length_chr;
            if start1 > chr_size; break; end
        end
        
        % Final trailing edge cleanup
        if start1 <= chr_size
            bins(bc,1) = c1;
            bins(bc,2) = start1;
            bins(bc,3) = chr_size;
            breaks_per_bin = find(sorted_pos(:,1) <= (chr_size + 1) & sorted_pos(:,1) >= start1);
            bins(bc,4) = length(breaks_per_bin);
            bc = bc + 1;
        end
    else
        % Fallback if chromosome is smaller than avg_bin_length
        if start1 <= chr_size
            bins(bc,1) = c1;
            bins(bc,2) = start1;
            bins(bc,3) = chr_size;
            breaks_per_bin = find(sorted_pos(:,1) <= (chr_size + 1) & sorted_pos(:,1) >= start1);
            bins(bc,4) = length(breaks_per_bin);
            bc = bc + 1;
        end
    end
end

numbins = size(bins,1);   
end