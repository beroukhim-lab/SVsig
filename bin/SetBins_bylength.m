function [bins, numbins] = SetBins_bylength(events,avg_bin_length,chsize,CHR,bin_shift)

% SetBins_bylength
% Build fixed-width genomic bins for the selected chromosomes.
%
% This function supports both the original unshifted lattice and translated
% lattices used for bin-shift runs. The important invariant is that the
% shift=0 case must reproduce the historical no-shift binning exactly.
%
% Inputs
%   events          Numeric event matrix. Columns 1/2 are chr/pos for the
%                   first breakpoint and columns 4/5 are chr/pos for the
%                   second breakpoint.
%   avg_bin_length  Target bin length used to derive a per-chromosome
%                   uniform lattice.
%   chsize          Vector of chromosome sizes indexed by chromosome id.
%   CHR             Chromosomes to include.
%   bin_shift       Optional non-negative shift in bp. When omitted or 0,
%                   bins are identical to the historical no-shift layout.
%
% Outputs
%   bins            N x 4 numeric matrix with columns:
%                     1 = chromosome
%                     2 = bin start (1-based, inclusive)
%                     3 = bin end   (1-based, inclusive)
%                     4 = breakpoint count in the bin
%   numbins         Number of rows in bins.

disp('setting bins boundaries...');

% Collect all breakpoint coordinates once so per-chromosome counting can be
% done from a sorted vector of positions.
chr = [events(:,1);events(:,4)];
pos = [events(:,2);events(:,5)];

bc = 1;
bins = [];

% Parse and validate the optional shift argument.
% The shift is defined in bp and is applied as a pure translation of the
% same base lattice used by the no-shift run.
if nargin >= 5 && ~isempty(bin_shift)
    shift = floor(double(bin_shift));
    if ~isfinite(shift) || shift < 0
        error('SetBins_bylength:BadBinShift', 'bin_shift must be a non-negative integer.');
    end
else
    shift = 0;
end

for c1 = CHR
    chr_size = floor(chsize(c1));
    if chr_size < 1
        continue;
    end

    % Restrict breakpoint counting to the current chromosome.
    chri = find(chr == c1);
    sorted_pos = sort(pos(chri));

    % Preserve the original no-shift geometry exactly.
    %
    % Historical behavior:
    %   num_chr        = floor(chr_size / avg_bin_length)
    %   avg_length_chr = floor(avg_bin_length + remainder/num_chr)
    %
    % That historical construction defines a uniform lattice with exactly
    % num_chr bins and no extra trailing cleanup bin. To keep shift=0 fully
    % compatible with legacy output, we keep that geometry unchanged and
    % implement shifted runs only as a translation of this same lattice.
    num_chr = floor(chr_size / avg_bin_length);
    if num_chr < 1
        % For very small chromosomes, fall back to a single bin spanning the
        % chromosome after clipping below.
        num_chr = 1;
    end
    avg_length_chr = floor(avg_bin_length + ((chr_size - num_chr * avg_bin_length) / num_chr));

    % Build the chromosome lattice one bin at a time.
    %
    % Unshifted lattice (shift = 0):
    %   bin k = [1 + (k-1)*L, k*L]
    %
    % Shifted lattice (shift > 0):
    %   bin k = [1 + shift + (k-1)*L, shift + k*L]
    %
    % where L = avg_length_chr.
    %
    % Importantly, shift=0 reduces algebraically to the original no-shift
    % bins. For shifted runs, we clip the translated bins to chromosome
    % bounds and skip bins that lie completely outside the chromosome.
    for c2 = 1:num_chr
        raw_start = 1 + shift + (c2 - 1) * avg_length_chr;
        raw_end   = shift + c2 * avg_length_chr;

        % Clip translated bin edges to the valid chromosome interval.
        start1 = max(1, raw_start);
        end1   = min(chr_size, raw_end);

        % If the translated bin falls completely beyond the chromosome, it
        % contributes nothing and should not be emitted.
        if start1 > chr_size || end1 < 1 || start1 > end1
            continue;
        end

        bins(bc,1) = c1;
        bins(bc,2) = start1;
        bins(bc,3) = end1;

        % Count breakpoints assigned to this bin.
        %
        % We intentionally preserve the original counting rule
        %   pos <= end + 1  and  pos >= start
        % to avoid changing downstream low-density filtering behavior.
        breaks_per_bin = find(sorted_pos(:,1) <= (end1 + 1) & sorted_pos(:,1) >= start1);
        bins(bc,4) = length(breaks_per_bin);

        bc = bc + 1;
    end
end

numbins = size(bins,1);
end
