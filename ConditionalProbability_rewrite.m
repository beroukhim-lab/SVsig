
function [sij,sij1dy] = ConditionalProbability_rewrite(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx)

disp ('calculating conditional probabilities...');
global firstbin 
global lastbin
global intra_chr_a


% setting up sij matrix
nume=length(events);
num_bins = length(bins);
sij = zeros(num_bins,num_bins,num_annot);
bpsize=sum(chsize(CHR));




for c1=CHR
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin;

%calculate sij matrix
for c2=firstbin:lastbin

       %
        sij(c2,c2+1:lastbin,:) = upper_diag;
        %calculate last bin
        last_diag = find(sij1dx<diag_bin_size,1,'last');
        
        
       %what exactly should conditional probability be?
       %prob of each tile / prob of tile that broke first?
       %then accounts for size of tile that is being invaded into? 
       %do we externally need to take into account some sort of length
       %factor? 
       
       %for length factor, take average distance between two bins? 
       
       %or if conditional probability, 
   end
    
    
    
    
    
    %inter chromsomal area in basepairs by chromosome * whole genome - area
    %of the inter chromosome
    inter_area = (lastbin-firstbin+1)*(bpsize-(bins(lastbin,3)-bins(firstbin,2)));
    for ca=1:num_annot
        sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca) + sij(chr_intra,chr_intra,ca)';
        %20201216 adding line to multiply intra by bin i/j
        sij(chr_intra,chr_intra,ca) = bsxfun(@times,sij(chr_intra,chr_intra,ca),(bins(chr_intra,3)-bins(chr_intra,2))');
       
        %normalize intrachromsomal events  
        intra_norm = intra_chr(c1,ca)*inter_area/(1-intra_chr(c1,ca))/sum(sum(sij(chr_intra,chr_intra,ca)));
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





