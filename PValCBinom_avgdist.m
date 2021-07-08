%takes in the binned matrix with weights and the background model
%probabilities 
%returns the qvalue that a tile is significant according to the
%beta_binomial distribution

function [qFDR, pa, pval_tophits,mfull] = PValCBinom_avgdist(mfull, p, qe, sij, bins, events)

global FDR_THRESHOLD

% general variables
mat_size = size(mfull);

% keep only upper diagonal part of mfull 
mfull = triu(full(mfull));
%divide diagonal elements by 2
mfull(eye(mat_size)~=0) = diag(mfull)/2;
%calculate N value (number of trials)  = total number of rearrangments
nume = sum(mfull(:));

% create probability matrix 
if isempty(p)
    pa = bsxfun(@times,qe,sij);
    pa = triu(pa + pa');
    pa(eye(mat_size)~=0) = diag(pa)/2;
else
    pa=p;
%if symmetric keep only the upper diagonal and divide diagonal elements by 2
    
    if issymmetric(p)
        pa=triu(p);
        pa(eye(mat_size)~=0) = diag(pa)/2;
    end
end

%double background rate as a filter
%pa = 2*pa; 

% divide tiles with non-zero values from tiles with zero values
high_k = find(mfull>=2 & pa>0);
pos_k = find(mfull==1 & pa>0);
zero_k = find(mfull==0 & pa>0);
mfull_pos = full(mfull(pos_k));
mfull_high = full(mfull(high_k));
p_high = pa(high_k);
p_pos = pa(pos_k);
p_zero = pa(zero_k);


%calculate p-values for each tile with >1 events according to continuous
%binomial 
disp(['calculating support for ' num2str(length(high_k)) ' tiles with >1 events']);
tic
t_dv=zeros(length(high_k),2);

for c1=1:length(high_k)
    [a1, a2]=ind2sub(mat_size,high_k(c1));     
    %[h, p_ks, n]=compare_avg_length( [a1 a2], bins, events, 0 );
    [p_mw, n]=compare_avg_length_complex([a1 a2], bins, events);

    t_dv(c1,1)=p_mw;
    t_dv(c1,2)=n;

end
%p_high_s=p_high;
toc
disp(['calculating p-val for ' num2str(length(high_k)) ' tiles with >1 events']);
tic
rand_nnz = rand(length(high_k),1);

%if ~approx_flag
%    pval_high = binopdf(mfull_high, nume, p_high).*rand_nnz+(1-binocdf(mfull_high,nume,p_high));
%else
    %pval_high = binopdf(mfull_high, nume, p_high).*rand_nnz+(1-binocdf(mfull_high,nume,p_high));
    %pval_high0 is pval_high without the random term
    %Note: pval_high0 > pval_high so the random term will push a hit into
    %significance
    %pval_high0 = (1-binocdf(mfull_high-1,nume,p_high));
    
pval_high0 = 1 - betacdf(p_high,mfull_high, nume + 1 - mfull_high,'upper');

%end
toc

%fisher's method to combine two pvals for pval_high0 terms 

combined_pval=zeros(length(high_k),1);
for c1=1:length(pval_high0)
    if t_dv(c1,2)>4
    chi_val = -2*(log(pval_high0(c1)) + log(t_dv(c1,1)));
    combined_pval(c1) = 1 - chi2cdf(chi_val,4); 
    else
    %combined_pval(c1) = pval_high0(c1); 
    chi_val = -2*(log(pval_high0(c1)) + log(t_dv(c1,1)));
    combined_pval(c1) = 1 - chi2cdf(chi_val,4);    
    end

end
    

%pvalue of zero event tiles is 1
pval_low = ones(1, length(zero_k))';
pval_pos=ones(1, length(pos_k))';



%pvalues for all the tiles
pval_tble=[pval_low zero_k;pval_pos pos_k;combined_pval high_k];
%pvalues for only the tiles with high counts without the random term
combined_pval=[combined_pval high_k];



pval=sortrows(pval_tble,1);
qFDR=mafdr(pval(:,1),'BHFDR','true');
%Kiran added FDR_THRESHOLD as a global variable
hits_idx=find(qFDR< FDR_THRESHOLD);
tophits=length(hits_idx);
pval_tophits = pval(1:tophits,:);
return 

