function [sij,sij1dy] = ConditionalProbability_mediandiag_interp(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx, bins_event_tble)

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
sij1dy(:,1,:) = EventLengthDist_G(sij1dx,events,0);
%line to account for higher sij1dx binning
%sij1dy(:,1,:)=(100/10).*sij1dy(:,1,:);
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

intra_events=events(events(:,1)==events(:,4),:);
median_interp_tble = zeros(length(intra_events),2);

%loop through all intra events (table of pairs of breakpoints)
for i = 1:length(intra_events)   
    %find the indices of the bins in which the ith and jth breaks fall in        
    bin1 = find(bins(:,1)==intra_events(i,1) & bins(:,2)<=intra_events(i,2) & bins(:,3)>=intra_events(i,2));
    bin2 = find(bins(:,1)==intra_events(i,4) & bins(:,2)<=intra_events(i,5) & bins(:,3)>=intra_events(i,5));
    
    median_interp_tble(i,1) = abs(((bins(bin1,3)+bins(bin1,2))/2) - ((bins(bin2,3)+bins(bin2,2))/2));
    median_interp_tble(i,2) = abs(intra_events(i,5) - intra_events(i,2));

end
%plot(log(sij1dx),log(sij1dy1000),'*k',log(sij1dx(1:99)),log(sij1dy1000interp(1:99)),'-k');
%plot(median_interp_tble(:,1),median_interp_tble(:,2),'*k');

[C,ia,idx] = unique(median_interp_tble(:,1),'stable');
val = accumarray(idx,median_interp_tble(:,2),[],@mean); 
unique_median_tble = [C val];

%plot(unique_median_tble(:,1),unique_median_tble(:,2),'*k');

           
%calculate sij matrix
for c1=CHR
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin; 
 
    
    for c2=firstbin:lastbin
        diag_bin = sum(bins(c2,2:3),2)/2;
        diag_bin_size=(bins(c2,3)-bins(c2,2));
        upper_diag_bins=sum(bins(c2:lastbin,2:3),2)/2-diag_bin;
        half_size=(bins(c2:lastbin,3)-bins(c2:lastbin,2))/2;
        length_interp=interp1(unique_median_tble(:,1),unique_median_tble(:,2),upper_diag_bins,'pchip');
        
        upper_diag=(interp1(sij1dx',squeeze(sij1dy(:,c1,:)),length_interp,'pchip'))/2;
        %sij(c2,c2+1:lastbin,:) = bsxfun(@times,upper_diag,bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2));
        sij(c2,c2:lastbin,:) = upper_diag;
        sij(c2,c2,:) = sij(c2,c2,:)/2;

        %sij(c2,c2+1:lastbin,:) = bsxfun(@times,upper_diag,bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2));

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





