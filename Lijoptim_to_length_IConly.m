function [sij1dy,opt1dy,fval,sij,intra_chr] = Lijoptim_to_length_IConly(events,chsize,bins,CHR,R0,mfull,sij1dx)
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
sij = ones(num_bins,num_bins);
bpsize=sum(chsize(CHR));

for c1=CHR,
    firstbin(c1)=find(bins(:,1)==c1,1);
    lastbin(c1)=find(bins(:,1)==c1,1,'last');
    chr_intra = firstbin(c1):lastbin(c1);
    b_mat(chr_intra,chr_intra)=1;
    
    sij(chr_intra,chr_intra) = zeros(length(chr_intra), length(chr_intra));
    
end


%intra_norm=(1-inter_chr)/sum(sum(Rij.*sij)); % here sij is still zero for all inter-chromosomal events
%sij=intra_norm*sij.*(b_mat==1)+inter_lij*(b_mat==0);
%sij = marginal_normalization( sij, R );

sij(find(eye(size(sij)))) = 0;


end







