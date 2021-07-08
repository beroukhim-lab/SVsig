function sij = sij_optimziation ( sij, p, qe, R )
%similar concept to marginal_normalization function for lijs 


%normalizes first row and column
norm_c=sij(1,:)*R;
sij(1,:)=sij(1,:)/norm_c;
sij(2:end,1)=sij(2:end,1)/norm_c;

%middle row/columns
for c1=2:length(sij)-1,
    norm_c=sij(c1,c1:end)*R(c1:end)/(1-sij(c1,1:c1-1)*R(1:c1-1));
    sij(c1,c1:end)=sij(c1,c1:end)/norm_c;
    sij(c1+1:end,c1)=sij(c1+1:end,c1)/norm_c;
end

%last tile
norm_c=sij(end,end)*R(end)/(1-sij(end,1:end-1)*R(1:end-1));
sij(end,end)=sij(end,end)/norm_c;

p0 = bsxfun(@times,qe,sij);
p=p0+p0'; 


qe = R/2; % qe=rand(length(R),1);qe=qe/sum(qe);
maxit=100;
%this essentially tries to solve A*x=b where A is sij + eye(size(sij)), b is R
%with default tolerance of 1e-6 (not specified here)
[qsolve,flag1,rr1,iter1] = pcg(sij+eye(size(sij)),R, [],maxit);



