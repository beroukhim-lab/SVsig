%takes in the binned matrix with weights and the background model
%probabilities 
%returns the qvalue that a tile is significant according to the
%beta_binomial distribution

function [qFDR, pa, pval_tophits,mfull] = PValCBinom(mfull, p, qe, sij)

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
pos_k = find(mfull>0 & pa>0);
zero_k = find(mfull==0 & pa>0);
mfull_pos = full(mfull(pos_k));
p_pos = pa(pos_k);
p_zero = pa(pos_k);


%calculate the p-values for each non-zero tile according to continous binomial and
%time the function

tic
disp('Calculating P Value using Continous Binomial')
%cdf_num = betainc(p_pos, mfull_pos, nume + 1 - mfull_pos, 'upper'); 
%cdf_denom = beta(nume + 1 - mfull_pos, mfull_pos);
%pval_nnz =  1 - cdf_num./cdf_denom;
%we can use the beta cdf to represent the continous binomial cdf
% 1 - cdfcbinom will give us 
pval_nnz = 1 - betacdf(p_pos,mfull_pos, nume + 1 - mfull_pos,'upper');
toc

%pvalue of zero event tiles is 1
pval_z = ones(1, length(zero_k))';
pval_tble=[pval_z zero_k;pval_nnz pos_k];


pval=sortrows(pval_tble,1);
qFDR=mafdr(pval(:,1),'BHFDR','true');
%Kiran added FDR_THRESHOLD as a global variable
hits_idx=find(qFDR< FDR_THRESHOLD);
tophits=length(hits_idx);
pval_tophits = pval(1:tophits,:);
return 

