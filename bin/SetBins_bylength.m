function [bins, numbins] = SetBins_bylength(events,avg_bin_length,chsize,CHR,bin_shift)

disp('setting bins boundaries...');

max_no_events=1e5;

chr = [events(:,1);events(:,4)];
pos = [events(:,2);events(:,5)];

bc=1;

use_shifted_bins = (nargin >= 5) && ~isempty(bin_shift);
if use_shifted_bins
    bin_shift = floor(double(bin_shift));
    if ~isfinite(bin_shift) || bin_shift < 1
        error('SetBins_bylength:BadBinShift', 'bin_shift must be a positive integer.');
    end
end

%for all the chrmosomes
for c1 = CHR

    chri = find(chr==c1);
    sorted_pos = sort(pos(chri));
    num_chr = floor(chsize(c1)/avg_bin_length);
    avg_length_chr = floor(avg_bin_length + ((chsize(c1)-num_chr*avg_bin_length)/num_chr));

    if ~use_shifted_bins
        % Legacy path kept to ensure n_binshifts=0 is identical to existing output.
        start1=1;
        for c2 = 1:num_chr
            bins(bc,1) = c1;
            bins(bc,2)=start1;
            bins(bc,3) =start1+avg_length_chr-1;
            breaks_per_bin=find(sorted_pos(:,1)<=(start1+avg_length_chr) & sorted_pos(:,1)>=start1);
            bins(bc,4) = length(breaks_per_bin);
            bc=bc+1;
            start1=start1+avg_length_chr;
        end
    else
        chr_size = floor(chsize(c1));

        % Shift bin at chromosome start: [0, bin_shift-1].
        shifted_end = min(bin_shift - 1, chr_size);
        bins(bc,1) = c1;
        bins(bc,2) = 0;
        bins(bc,3) = shifted_end;
        breaks_per_bin = find(sorted_pos(:,1) <= shifted_end & sorted_pos(:,1) >= 0);
        bins(bc,4) = length(breaks_per_bin);
        bc = bc + 1;

        start1 = bin_shift;
        while start1 <= chr_size
            end1 = min(start1 + avg_length_chr - 1, chr_size);
            bins(bc,1) = c1;
            bins(bc,2) = start1;
            bins(bc,3) = end1;
            breaks_per_bin = find(sorted_pos(:,1) <= end1 & sorted_pos(:,1) >= start1);
            bins(bc,4) = length(breaks_per_bin);
            bc = bc + 1;
            start1 = start1 + avg_length_chr;
        end
    end

end

numbins = size(bins,1);   

