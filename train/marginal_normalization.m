function sij = marginal_normalization( sij, R )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
