function [sij1dy,opt1dy,fval,sij,intra_chr] = Lijoptim_to_length_test(events,chsize,bins,CHR,R0,mfull,sij1dx)
% finds maximum likelehood 1D distribution for mult-model 

%len_factor = model_length_dist( events,bins, CHR, sij1dx);

intra_chr=zeros(length(CHR),length(CHR));
%R=R0.*(bins(:,3)-bins(:,2));
R=R0/2;
R=R/sum(R);
Rij=kron(R,R');
mat_size=size(Rij);

b_mat=zeros(mat_size);
nume=length(events);

for c1=CHR,
    bin_ind(:,c1)=bins(:,1)==c1;
end
inter_chr=0;
inter_rij=0;
for c1=CHR,
    for c2=CHR,
        if c1~=c2
            inter_chr=inter_chr+sum(sum(mfull(bin_ind(:,c1),bin_ind(:,c2))));
            inter_rij=inter_rij+sum(sum(Rij(bin_ind(:,c1),bin_ind(:,c2))));
        end
    end
end
inter_chr=inter_chr/2/nume;
inter_lij=inter_chr/inter_rij;

% calculate the 1D Sij distribution
%sidjx indicates bins that are chosen for empirical distribution
sij1dy = EventLengthDist_G(sij1dx,events,0);

%to test changes in 1/L distribution
n=1;
sij1dy=(sij1dy).^n;

%to account for differing sij1dy values
%sij1dy=(100/10).*sij1dy;

%sij1dy = sum(sij1dy,2)/length(sij1dy(1,:));     %averaging across types of events (4)
%sij1dy = sij1dy/sum(sij1dy);                    %normalizing to sum to 1


%fixing sij1dy average to account for event type weight
nume=length(events)
annot_frac(1) = sum(events(:,3)==1&events(:,6)==1)/nume;
annot_frac(2) = sum(events(:,3)==1&events(:,6)==2)/nume;
annot_frac(3) = sum(events(:,3)==2&events(:,6)==1)/nume;
annot_frac(4) = sum(events(:,3)==2&events(:,6)==2)/nume;
for ca=1:4
    sij1dy(ca,:) = sij1dy(ca,:) * annot_frac(ca);
end
sij1dy = sum(sij1dy,2);     % summing across types of events (4)
sij1dy = sij1dy/sum(sij1dy);           %normalizing to sum to 1

log_sij1dy=log10(sij1dy(1:end-1));
sij1dx=sij1dx';
d_sij1dx=diff(sij1dx);
opt1dy=[];
fval=[];

% calculate sij 
num_bins = length(bins);
sij = zeros(num_bins,num_bins);
bpsize=sum(chsize(CHR));

for c1=CHR,
    firstbin(c1)=find(bins(:,1)==c1,1);
    lastbin(c1)=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin(c1):lastbin(c1);
    b_mat(chr_intra,chr_intra)=1;
    
    for c2=firstbin(c1):lastbin(c1),
        diag_bin(c2) = sum(bins(c2,2:3),2)/2;
        %diag_bin_size is equivalent to half_size in ConditionalProbability.m
        diag_bin_size(c2)=(bins(c2,3)-bins(c2,2))/2;
        %does cubic interpolation for start and end of each bin (based on bp
        %distance away from diag bin), averages these starts and ends 
        %seems to be essentially the same as for sijs, written differently
        sij(c2,c2+1:lastbin(c1)) = (interp1(sij1dx,sij1dy,bins(c2+1:lastbin(c1),2)-diag_bin(c2),'pchip')+interp1(sij1dx,sij1dy,bins(c2+1:lastbin(c1),3)-diag_bin(c2),'pchip'))'/2;
        last_diag(c2) = find(sij1dx<diag_bin_size(c2),1,'last');
        sij(c2,c2) = ( sij1dy(1:last_diag(c2)-1)'*d_sij1dx(1:last_diag(c2)-1) + (diag_bin_size(c2) - sij1dx(last_diag(c2))-1) * interp1(sij1dx,sij1dy,diag_bin_size(c2),'pchip'))/(diag_bin_size(c2)-sij1dx(1));  
        %sij(c2,c2) = ( sij1dy(1:last_diag(c2)-1)'*d_sij1dx(1:last_diag(c2)-1) + (diag_bin_size(c2) - sij1dx(last_diag(c2))-1) * interp1(sij1dx,sij1dy,diag_bin_size(c2),'pchip'))/(diag_bin_size(c2));  
        sij(c2, c2) = 1;

    end
    sij(chr_intra,chr_intra) = sij(chr_intra,chr_intra) + sij(chr_intra,chr_intra)';
end


intra_norm=(1-inter_chr)/sum(sum(Rij.*sij)); % here sij is still zero for all inter-chromosomal events
sij=intra_norm*sij.*(b_mat==1)+inter_lij*(b_mat==0);

%add here to normalize lij by event ratios before marg norm lol this is so janky
[sij]=renormalize_tiles(mfull, sij, events, bins, CHR);

sij = marginal_normalization( sij, R );

end








