function [sij,sij1dy] = ConditionalProbability_copy(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx)

disp('calculating conditional probabilities...');
global firstbin
global lastbin
global intra_chr_a

nChr = numel(CHR);

% calculate the intra-chromosomal event fraction from data
intra_chr = zeros(nChr,num_annot);

if 0 % intra_chr ratio per chromosome
    for k = 1:nChr
        cChr = CHR(k);
        bin_ind = bins(:,1)==cChr;

        for c2 = 1:num_annot
            denom = sum(sum(mfull{c2}(bin_ind,:)));
            if denom == 0
                intra_chr(k,c2) = 0;
            else
                intra_chr(k,c2) = sum(sum(mfull{c2}(bin_ind,bin_ind))) / denom;
            end
        end
    end
else % global intra_chr ratio
    intra_chr_a = zeros(1,num_annot);

    for c2 = 1:num_annot
        for k = 1:nChr
            cChr = CHR(k);
            bin_ind = bins(:,1)==cChr;
            intra_chr_a(c2) = intra_chr_a(c2) + sum(sum(mfull{c2}(bin_ind,bin_ind)));
        end

        denom = sum(sum(mfull{c2}(:,:)));
        if denom == 0
            intra_chr_a(c2) = 0;
        else
            intra_chr_a(c2) = intra_chr_a(c2) / denom;
        end
    end

    intra_chr = repmat(intra_chr_a,nChr,1);
end

nume = size(events,1);

% calculate the 1D Sij distribution
sij1dy = zeros(length(sij1dx),1,num_annot);
sij1dy(:,1,:) = EventLengthDist_G(sij1dx,events,0);

% repeat genome-wide distribution for each chromosome included in CHR
sij1dy = repmat(sij1dy,1,nChr,1);

% calculate sij
num_bins = size(bins,1);
sij = zeros(num_bins,num_bins,num_annot);
bpsize = sum(chsize(CHR));

% distances between sij1dx bins
d_sij1dx = diff(sij1dx)';

sd_sij1dx = zeros(1,length(d_sij1dx));
sd_sij1dx(1) = d_sij1dx(1)/2;

for c1 = 2:length(d_sij1dx)
    sd_sij1dx(1,c1) = sd_sij1dx(c1-1) + sum(d_sij1dx(c1-1:c1))/2;
end

% IMPORTANT:
% sij1dy second dimension is 1:nChr, not chromosome labels.
% So index it with k, not cChr.
sij1_area = zeros(length(sij1dx)-1,nChr,num_annot);

for k = 1:nChr
    sij1_area(:,k,:) = bsxfun(@times, squeeze(sij1dy(1:end-1,k,:)), d_sij1dx);
end

% calculate sij matrix
for k = 1:nChr

    cChr = CHR(k);

    firstbin = find(bins(:,1)==cChr,1);
    lastbin  = find(bins(:,1)==cChr,1,'last');

    if isempty(firstbin) || isempty(lastbin)
        warning('No bins found for chromosome %g. Skipping.', cChr);
        continue;
    end

    chr_intra = firstbin:lastbin;

    % If only one bin in chromosome, avoid c2+1:lastbin empty-edge issues
    if firstbin == lastbin
        c2 = firstbin;
        sij(c2,c2,:) = 1;

        for ca = 1:num_annot
            sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca) + ...
                                          sij(chr_intra,chr_intra,ca)';

            if lastbin < num_bins
                cols_right = lastbin+1:num_bins;
                sij(firstbin:lastbin,cols_right,ca) = ...
                    repmat((bins(cols_right,3)-bins(cols_right,2))', ...
                    lastbin-firstbin+1,1);
            end

            if firstbin > 1
                cols_left = 1:firstbin-1;
                sij(firstbin:lastbin,cols_left,ca) = ...
                    repmat((bins(cols_left,3)-bins(cols_left,2))', ...
                    lastbin-firstbin+1,1);
            end
        end

        continue;
    end

    % First bin in chromosome
    c2 = firstbin;

    diag_bin = sum(bins(c2,2:3),2)/2;
    diag_bin_size = bins(c2,3)-bins(c2,2);

    upper_diag_bins = sum(bins(c2+1:lastbin,2:3),2)/2 - diag_bin;
    half_size = (bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2))/2;

    upper_diag = ( ...
        interp1(sij1dx', squeeze(sij1dy(:,k,:)), upper_diag_bins-half_size, 'pchip') + ...
        interp1(sij1dx', squeeze(sij1dy(:,k,:)), upper_diag_bins+half_size, 'pchip') ...
        ) / 2;

    sij(c2,c2+1:lastbin,:) = upper_diag;

    last_diag = find(sij1dx < diag_bin_size,1,'last');

    if isempty(last_diag) || last_diag <= 1
        sij(c2,c2,:) = 1;
    elseif last_diag == 2
        sij(c2,c2,:) = ( ...
            (1 - sd_sij1dx(1:last_diag-1)/diag_bin_size) * ...
            squeeze(sij1_area(1:last_diag-1,k,:))' + ...
            (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * ...
            interp1(sij1dx', squeeze(sij1dy(:,k,:)), diag_bin_size, 'pchip') ...
            );
    else
        sij(c2,c2,:) = ( ...
            (1 - sd_sij1dx(1:last_diag-1)/diag_bin_size) * ...
            squeeze(sij1_area(1:last_diag-1,k,:)) + ...
            (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * ...
            interp1(sij1dx', squeeze(sij1dy(:,k,:)), diag_bin_size, 'pchip') ...
            );
    end

    sij(c2,c2,:) = 1;

    % Remaining bins in chromosome
    for c2 = firstbin+1:lastbin

        diag_bin = sum(bins(c2,2:3),2)/2;
        diag_bin_size = bins(c2,3)-bins(c2,2);

        if c2 < lastbin
            upper_diag_bins = sum(bins(c2+1:lastbin,2:3),2)/2 - diag_bin;
            half_size = (bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2))/2;

            upper_diag = ( ...
                interp1(sij1dx', squeeze(sij1dy(:,k,:)), upper_diag_bins-half_size, 'pchip') + ...
                interp1(sij1dx', squeeze(sij1dy(:,k,:)), upper_diag_bins+half_size, 'pchip') ...
                ) / 2;

            sij(c2,c2+1:lastbin,:) = upper_diag;
        end

        last_diag = find(sij1dx < diag_bin_size,1,'last');

        if isempty(last_diag) || last_diag <= 1
            sij(c2,c2,:) = 1;
        elseif last_diag == 2
            sij(c2,c2,:) = ( ...
                (1 - sd_sij1dx(1:last_diag-1)/diag_bin_size) * ...
                squeeze(sij1_area(1:last_diag-1,k,:))' + ...
                (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * ...
                interp1(sij1dx', squeeze(sij1dy(:,k,:)), diag_bin_size, 'pchip') ...
                );
        else
            sij(c2,c2,:) = ( ...
                (1 - sd_sij1dx(1:last_diag-1)/diag_bin_size) * ...
                squeeze(sij1_area(1:last_diag-1,k,:)) + ...
                (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * ...
                interp1(sij1dx', squeeze(sij1dy(:,k,:)), diag_bin_size, 'pchip') ...
                );
        end

        sij(c2,c2,:) = bsxfun(@times, sij(c2,c2,:), 1/(bins(c2,3)-bins(c2,2)));
        sij(c2,c2,:) = 1;
    end

    % Inter-chromosomal weights
    inter_area = (lastbin-firstbin+1) * ...
                 (bpsize - (bins(lastbin,3)-bins(firstbin,2)));

    for ca = 1:num_annot

        sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca) + ...
                                      sij(chr_intra,chr_intra,ca)';

        if lastbin < num_bins
            sij(firstbin:lastbin,lastbin+1:end,ca) = ...
                repmat((bins(lastbin+1:end,3)-bins(lastbin+1:end,2))', ...
                lastbin-firstbin+1,1);
        end

        if firstbin > 1
            sij(firstbin:lastbin,1:firstbin-1,ca) = ...
                repmat((bins(1:firstbin-1,3)-bins(1:firstbin-1,2))', ...
                lastbin-firstbin+1,1);
        end
    end
end

annot_frac(1) = sum(events(:,3)==1 & events(:,6)==1) / nume;
annot_frac(2) = sum(events(:,3)==1 & events(:,6)==2) / nume;
annot_frac(3) = sum(events(:,3)==2 & events(:,6)==1) / nume;
annot_frac(4) = sum(events(:,3)==2 & events(:,6)==2) / nume;

for ca = 1:num_annot

    sij(:,:,ca) = renormalize_tiles(mfull{ca}, sij(:,:,ca), events, bins, CHR);

    row_sums = sum(sij(:,:,ca),2);
    row_sums(row_sums == 0) = 1;

    sij(:,:,ca) = bsxfun(@rdivide, sij(:,:,ca), row_sums);
    sij(:,:,ca) = sij(:,:,ca) * annot_frac(ca);

end

% average intra_chr across all rearrangement types
intra_chr_a = mean(intra_chr_a);

end
