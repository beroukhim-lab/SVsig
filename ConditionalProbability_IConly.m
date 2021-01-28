function [sij,sij1dy] = ConditionalProbability(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx)

disp ('calculating conditional probabilities...');
global firstbin 
global lastbin
global intra_chr_a


%calculate the intra-chromosomal event fraction from data
intra_chr=zeros(length(CHR),num_annot);
if 0 % intra_chr ratio pre chromosome
    for c1=CHR,
        bin_ind=bins(:,1)==c1;
        for c2=1:num_annot,
            intra_chr(c1,c2)=sum(sum(mfull{c2}(bin_ind,bin_ind)))/sum(sum(mfull{c2}(bin_ind,:)));
        end
    end
else  % global intra_chr ratio 
    intra_chr_a=zeros(1,num_annot);
    for c2=1:num_annot,
        for c1=CHR,
            bin_ind=bins(:,1)==c1;
            intra_chr_a(c2)=intra_chr_a(c2)+sum(sum(mfull{c2}(bin_ind,bin_ind)));
        end
        intra_chr_a(c2)=intra_chr_a(c2)/sum(sum(mfull{c2}(:,:)));
    end
    intra_chr=repmat(intra_chr_a,length(CHR),1);
end
nume=length(events);

% calculate the 1D Sij distribution
% sij1dy = EventLengthDist(sij1dx,events,chsize,EventLengthThreshold,CHR,0);
% sij1dy0=sum(sij1dy,2)/length(sij1dy(1,:,1));
% for c1=1:length(sij1dy(1,:,1)),
%     sij1dy(:,c1,:)=sij1dy0;
% end

sij1dy = zeros(length(sij1dx),1,num_annot);
%UNCOMMENT HERE
%sij1dy(:,1,:) = EventLengthDist_G(sij1dx,events,0);
%returns genome-wide avergae number of events per bp by length for each
%SV-type, grouped by sijdx boundaries (repeated for all chr)
sij1dy=repmat(sij1dy,1,length(CHR),1);

% calculate sij 
num_bins = length(bins);
sij = zeros(num_bins,num_bins,num_annot);
bpsize=sum(chsize(CHR));


%sd_sij1dx seems to be average of adjacent d_sij1dx values? but why? 
d_sij1dx=diff(sij1dx)';
sd_sij1dx(1)=d_sij1dx(1)/2;
for c1=2:length(d_sij1dx)
    sd_sij1dx(1,c1)=sd_sij1dx(c1-1)+sum(d_sij1dx(c1-1:c1))/2;
end

%UNCOMMENT NEXT 3 ROWS
%for c1=CHR
%    sij1_area(:,c1,:)=(bsxfun(@times,squeeze(sij1dy(1:end-1,c1,:)),d_sij1dx));
%end


%calculate sij matrix
for c1=CHR
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin; 
      
    %inter chromsomal area in basepairs by chromosome * whole genome - area
    %of the inter chromosome
    inter_area = (lastbin-firstbin+1)*(bpsize-(bins(lastbin,3)-bins(firstbin,2)));
    for ca=1:num_annot
        sij(chr_intra,chr_intra,ca) = repmat(0, length(chr_intra), length(chr_intra));
        %20201216 adding line to multiply intra by bin i/j
        %sij(chr_intra,chr_intra,ca) = bsxfun(@times,sij(chr_intra,chr_intra,ca),(bins(chr_intra,3)-bins(chr_intra,2))');
        %normalize intrachromsomal events  
        intra_norm = 0;
        sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca)*intra_norm;
        %giving value of 1 to all inter_chromosomal bins
        sij(firstbin:lastbin,lastbin+1:end,ca) = repmat((bins(lastbin+1:end,3)-bins(lastbin+1:end,2))',lastbin-firstbin+1,1);
        sij(firstbin:lastbin,1:firstbin-1,ca) = repmat((bins(1:firstbin-1,3)-bins(1:firstbin-1,2))',lastbin-firstbin+1,1);
    end
end

annot_frac(1) = sum(events(:,3)==1&events(:,6)==1)/nume;
annot_frac(2) = sum(events(:,3)==1&events(:,6)==2)/nume;
annot_frac(3) = sum(events(:,3)==2&events(:,6)==1)/nume;
annot_frac(4) = sum(events(:,3)==2&events(:,6)==2)/nume;
for ca=1:num_annot
    sij(:,:,ca) = bsxfun(@rdivide,sij(:,:,ca),sum(sij(:,:,ca), 2));
    sij(:,:,ca) = sij(:,:,ca) * annot_frac(ca);
end

%average intra_chr across all 4 rearrangement types
intra_chr_a = mean(intra_chr_a);





