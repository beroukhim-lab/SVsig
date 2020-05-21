function [sij1dy,opt1dy,fval,sij,intra_chr] = Lijoptim_to_length(events,chsize,bins,CHR,R0,mfull,sij1dx,len_factor)
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

sij1dy = EventLengthDist_G(sij1dx,events,0);
sij1dy = sum(sij1dy,2)/length(sij1dy(1,:));
sij1dy = sij1dy/sum(sij1dy);
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
        diag_bin_size(c2)=(bins(c2,3)-bins(c2,2))/2;
        sij(c2,c2+1:lastbin(c1)) = (interp1(sij1dx,sij1dy,bins(c2+1:lastbin(c1),2)-diag_bin(c2),'pchip')+interp1(sij1dx,sij1dy,bins(c2+1:lastbin(c1),3)-diag_bin(c2),'pchip'))'/2;
        last_diag(c2) = find(sij1dx<diag_bin_size(c2),1,'last');
        sij(c2,c2) = ( sij1dy(1:last_diag(c2)-1)'*d_sij1dx(1:last_diag(c2)-1) + (diag_bin_size(c2) - sij1dx(last_diag(c2))-1) * interp1(sij1dx,sij1dy,diag_bin_size(c2),'pchip'))/(diag_bin_size(c2)-sij1dx(1));    
    end
    sij(chr_intra,chr_intra) = sij(chr_intra,chr_intra) + sij(chr_intra,chr_intra)';
end

intra_norm=(1-inter_chr)/sum(sum(Rij.*sij)); % here sij is still zero for all inter-chromosomal events
sij=intra_norm*sij.*(b_mat==1)+inter_lij*(b_mat==0);
sij = marginal_normalization( sij, R );

%comment out after here to do without normalization
mfull=triu(mfull);
mfull(eye(mat_size)==1)=diag(mfull)/2;


len1d=length(sij1dy)-1;
ct=0;

sij1dy_init=log10(sij1dy(1:end-1));




%sij1dy_init=-10*rand(len1d,1);

options = optimoptions('fmincon','DiffMinChange',1e-3,'TolFun',1e-1,'TolX',1e-4);
% %[opt1dy, fval] = fmincon(@lij_optim_fun,sij1dy_init,[],[],repmat(diff(sij1dx)',len1d,1),ones(len1d,1),zeros(len1d,1),ones(len1d,1),[],options);
[opt1dy, fval] = fmincon(@lij_optim_fun,sij1dy_init,[],[],[],[],-20*ones(len1d,1),zeros(len1d,1),[],options);
% %[opt1dy, fval] = fminunc(@lij_optim_fun,sij1dy_init,options);
opt1dy = [10.^opt1dy;0];
opt1dy(1:end-1)=opt1dy(1:end-1)./sum(opt1dy(1:end-1).*diff(sij1dx))/nume*(sum(mfull(:)));

    function err=lij_optim_fun(log_sij1d)

        sij1d=[10.^(log_sij1d);0];
% %        sij1d=[log_sij1d;0];
        sij = zeros(num_bins,num_bins);
        for c1=CHR,
            chr_intra = firstbin(c1):lastbin(c1);
            for c2=firstbin(c1):lastbin(c1),
                sij(c2,c2+1:lastbin(c1)) = (interp1(sij1dx,sij1d,bins(c2+1:lastbin(c1),2)-diag_bin(c2),'pchip')+interp1(sij1dx,sij1d,bins(c2+1:lastbin(c1),3)-diag_bin(c2),'pchip'))'/2;
                sij(c2,c2) = ( sij1d(1:last_diag(c2)-1)'*d_sij1dx(1:last_diag(c2)-1) + (diag_bin_size(c2) - sij1dx(last_diag(c2))-1) * interp1(sij1dx,sij1d,diag_bin_size(c2),'pchip'))/(diag_bin_size(c2)-sij1dx(1));          
            end
            sij(chr_intra,chr_intra) = sij(chr_intra,chr_intra) + sij(chr_intra,chr_intra)';
        end

        intra_norm=(1-inter_chr)/sum(sum(Rij.*sij)); % here sij is still zero for all inter-chromosomal events
        sij=intra_norm*sij.*(b_mat==1)+inter_lij*(b_mat==0);
        sij = marginal_normalization( sij, R );
        sij = sij.*(sij>0);

        p_mult = Rij.*sij;
        p_mult = triu(2*p_mult/sum(p_mult(:)));
        p_mult(eye(mat_size)~=0)=diag(p_mult)/2;
        
        p_mult=p_mult';
        p_vec=p_mult(logical(tril(b_mat)));

        sij1dy_model=zeros(length(sij1dy)-1);
        for c1=1:length(len_factor(1,:)),
            sij1dy_model(c1)=p_vec'*len_factor(:,c1)/d_sij1dx(c1);
        end
        sij1dy_model=sij1dy_model/sum(sij1dy_model);
        
        err=sum(abs(log10(sij1dy_model)-log_sij1dy).^2);
        
        if mod(ct,len1d)==0;
            err
        end
        ct=ct+1;
%         %          figure(1)
% %          loglog(sij1dx,abs(sij1d-sij1dy),'o')
    end

end







