function sij = ConditionalProbability_fragile(events,fragile_event_annot,bins_event_tble,chsize,bins,CHR,sij1dx)

disp ('calculating conditional probabilities...');

% general vars
nume = length(events);
num_bins = length(bins);
sij = zeros(num_bins,num_bins);
bpsize=sum(chsize(CHR));
d_sij1dx=diff(sij1dx)';

%calculate the intra-chromosomal event fraction from data
intra_chr=sum(events(:,1)==events(:,4))/nume;

% get gloabl length distribution from data
sij1dy_g = EventLengthDist_G(sij1dx,events(~fragile_event_annot,:),0);
sij1dy_g = sum(sij1dy_g,2);
sij1dy_g(end-1)=sij1dy_g(end-1)*3;
sij1dy_g = sij1dy_g./sum(sij1dy_g(1:end-1).*d_sij1dx);
sij1dycdf_g=sij1dy_g(1:end-1).*d_sij1dx;

% get length distribution for fragile sites
sij1dy_f = EventLengthDist_G(sij1dx,events(logical(fragile_event_annot),:),0);
sij1dy_f = sum(sij1dy_f,2);
sij1dy_f = sij1dy_f./sum(sij1dy_f(1:end-1).*d_sij1dx);
sij1dycdf_f=sij1dy_f(1:end-1).*d_sij1dx;

% find bins with majority of fragile sites
fragile_bins=annotate_bins(bins,bins_event_tble,fragile_event_annot,0);

% calcualte sij
sd_sij1dx(1)=d_sij1dx(1)/2;
for c1=2:length(d_sij1dx),
    sd_sij1dx(1,c1)=sd_sij1dx(c1-1)+sum(d_sij1dx(c1-1:c1))/2;
end


for c1=CHR,
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin; 
    
    for c2=firstbin:lastbin,
        if fragile_bins(c2)
            sij1dy=sij1dy_f;
            sij1dycdf=sij1dycdf_f;
        else
            sij1dy=sij1dy_g;
            sij1dycdf=sij1dycdf_g;
        end
        diag_bin = sum(bins(c2,2:3),2)/2;
        diag_bin_size=(bins(c2,3)-bins(c2,2));
        diag_bins=abs(sum(bins(firstbin:lastbin,2:3),2)/2-diag_bin);
        half_size=(bins(firstbin:lastbin,3)-bins(firstbin:lastbin,2))/2;
        upper_diag=(interp1(sij1dx',sij1dy,diag_bins-half_size,'pchip')+interp1(sij1dx',sij1dy,diag_bins+half_size,'pchip'))/2;
%        upper_diag=(10.^interp1(log10(sij1dx'),log10(sij1dy),log10(abs(diag_bins-half_size)),'pchip')+10.^interp1(log10(sij1dx'),log10(sij1dy),log10(diag_bins+half_size),'pchip'))/2;
        sij(c2,firstbin:lastbin) = upper_diag.*(bins(firstbin:lastbin,3)-bins(firstbin:lastbin,2));
        last_diag = find(sij1dx<diag_bin_size,1,'last');
        sij(c2,c2) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*sij1dycdf(1:last_diag-1) + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',sij1dy,diag_bin_size,'pchip'));
    end
    
    % normalization factor set by inter- to intra- number of events ratio
    inter_area = (lastbin-firstbin+1)*(bpsize-(bins(lastbin,3)-bins(firstbin,2)));
    sij(chr_intra,chr_intra) = sij(chr_intra,chr_intra) + sij(chr_intra,chr_intra)';
    intra_norm = intra_chr*inter_area/(1-intra_chr)/sum(sum(sij(chr_intra,chr_intra)));
    sij(chr_intra,chr_intra) = sij(chr_intra,chr_intra)*intra_norm;
    sij(firstbin:lastbin,lastbin+1:end) = repmat((bins(lastbin+1:end,3)-bins(lastbin+1:end,2))',lastbin-firstbin+1,1);
    sij(firstbin:lastbin,1:firstbin-1) = repmat((bins(1:firstbin-1,3)-bins(1:firstbin-1,2))',lastbin-firstbin+1,1);
end

sij(:,:) = bsxfun(@rdivide,sij,sum(sij,2));

end





