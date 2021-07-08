function [sij,sij1dy] = ConditionalProbability_formula(events,chsize,bins,EventLengthThreshold,CHR,num_annot,mfull,sij1dx)

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
for c1=CHR
    sij1_area(:,c1,:)=(bsxfun(@times,squeeze(sij1dy(1:end-1,c1,:)),d_sij1dx));
end

count1=0;
count2=0;
count3=0;
count4=0;

%calculate sij matrix
for c1=CHR
    firstbin=find(bins(:,1)==c1,1);
    lastbin=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin:lastbin; 
    
    c2=firstbin;
    %find midpoint of bin
    diag_bin = sum(bins(c2,2:3),2)/2;
    diag_bin_size=(bins(c2,3)-bins(c2,2));
    %find distance of other bins in row from diag bin
    upper_diag_bins=sum(bins(c2+1:lastbin,2:3),2)/2-diag_bin;
    %bp size of each bin in the row divided by 2 
    half_size=(bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2))/2;
    %does cubic interpolation for start and end of each bin (based on bp
    %distance away from diag bin), averages these starts and ends
    upper_diag=(interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins-half_size,'pchip')+interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins+half_size,'pchip'))/2;
    %sij(c2,c2+1:lastbin,:) = bsxfun(@times,upper_diag,bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2));
    sij(c2,c2+1:lastbin,:) = upper_diag;
    last_diag = find(sij1dx<diag_bin_size,1,'last');
    % add similar correction to bin 2 and 3...
    %Kiran: fix bug for last_diag = = 2 
        if last_diag == 2 
        sij(c2, c2, :) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*squeeze(sij1_area(1:last_diag-1,c1,:))' + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));      
        else
        disp(c1)
        sij(c2,c2,:) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*squeeze(sij1_area(1:last_diag-1,c1,:)) + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));
        end 
        sij(c2, c2, :) = bsxfun(@times,sij(c2, c2, :), 1/(bins(c2,3)-bins(c2,2)));
        count1=count1+ sij(c2, c2, 1) ;
        count2=count2+ sij(c2, c2, 2) ;
        count3=count3+ sij(c2, c2, 3) ;
        count4=count4+ sij(c2, c2, 4) ;

   % sij(c2,c2,:) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*squeeze(sij1_area(1:last_diag-1,c1,:)) + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));
%    sij(c2,c2,:) = (sum(bsxfun(@times,squeeze(sij1dy(1:last_diag-1,c1,:)),diff(sij1dx(1:(last_diag)))')) + (diag_bin_size - sij1dx(last_diag)-1) * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));
    for c2=firstbin+1:lastbin
        diag_bin = sum(bins(c2,2:3),2)/2;
        diag_bin_size=(bins(c2,3)-bins(c2,2));
        upper_diag_bins=sum(bins(c2+1:lastbin,2:3),2)/2-diag_bin;
        half_size=(bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2))/2;
        upper_diag=(interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins-half_size,'pchip')+interp1(sij1dx',squeeze(sij1dy(:,c1,:)),upper_diag_bins+half_size,'pchip'))/2;
        %sij(c2,c2+1:lastbin,:) = bsxfun(@times,upper_diag,bins(c2+1:lastbin,3)-bins(c2+1:lastbin,2));
        sij(c2,c2+1:lastbin,:) = upper_diag;
        last_diag = find(sij1dx<diag_bin_size,1,'last');
        %Kiran: fix bug for last_diag = = 2 
        if last_diag == 2 
        sij(c2, c2, :) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*squeeze(sij1_area(1:last_diag-1,c1,:))' + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));      
        else
        sij(c2,c2,:) = ( (1-sd_sij1dx(1:last_diag-1)/diag_bin_size)*squeeze(sij1_area(1:last_diag-1,c1,:)) + (diag_bin_size - sij1dx(last_diag)-1)^2/diag_bin_size/2 * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));
        end  
        sij(c2, c2, :) = bsxfun(@times,sij(c2, c2, :), 1/(bins(c2,3)-bins(c2,2)));
        
        count1=count1+ sij(c2, c2, 1) ;
        count2=count2+ sij(c2, c2, 2) ;
        count3=count3+ sij(c2, c2, 3) ;
        count4=count4+ sij(c2, c2, 4) ;


        %        sij(c2,c2,:) = (sum(bsxfun(@times,squeeze(sij1dy(1:last_diag-1,c1,:)),diff(sij1dx(1:(last_diag)))')) + (diag_bin_size - sij1dx(last_diag)-1) * interp1(sij1dx',squeeze(sij1dy(:,c1,:)),diag_bin_size,'pchip'));
    end
    

    
    %inter chromsomal area in basepairs by chromosome * whole genome - area
    %of the inter chromosome
    inter_area = (lastbin-firstbin+1)*(bpsize-(bins(lastbin,3)-bins(firstbin,2)));
    for ca=1:num_annot
        sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca) + sij(chr_intra,chr_intra,ca)';
        %20201216 adding line to multiply intra by bin i/j
        %sij(chr_intra,chr_intra,ca) = bsxfun(@times,sij(chr_intra,chr_intra,ca),(bins(chr_intra,3)-bins(chr_intra,2))');

        intra1=mean(intra_chr_a);
        intrasum=sum(sij(chr_intra,chr_intra,ca), 'all');
        interavg=(intrasum)*((1-intra1)/intra1);
        tiles1=numel(sij(firstbin:lastbin,lastbin+1:end,ca));
        tiles2=numel(sij(firstbin:lastbin,1:firstbin-1,ca));

        interavg=interavg/(tiles1+tiles2);

        %normalize intrachromsomal events  
        %intra_norm = intra_chr(c1,ca)*inter_area/(1-intra_chr(c1,ca))/sum(sum(sij(chr_intra,chr_intra,ca)));
        %sij(chr_intra,chr_intra,ca) = sij(chr_intra,chr_intra,ca)*intra_norm;
        %giving value of 1 to all inter_chromosomal bins
        %sij(firstbin:lastbin,lastbin+1:end,ca) = repmat((bins(lastbin+1:end,3)-bins(lastbin+1:end,2))',lastbin-firstbin+1,1);
        %sij(firstbin:lastbin,1:firstbin-1,ca) = repmat((bins(1:firstbin-1,3)-bins(1:firstbin-1,2))',lastbin-firstbin+1,1);
        sij(firstbin:lastbin,lastbin+1:end,ca) = repmat(interavg, size(sij(firstbin:lastbin,lastbin+1:end,ca)));
        sij(firstbin:lastbin,1:firstbin-1,ca) = repmat(interavg, size(sij(firstbin:lastbin,1:firstbin-1,ca)));
       
        
 
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





