function [mix_model,opt_alpha, f_bic] = mix_model_param1( mfull, model1, model2, events, bins, CHR)
global complex
global weights


%model 1 is break invasion
%model 2 is double break join
%alpha represents the fraction of break invasion
% finds the alpha's of the mix model pij = alpha*model1 + (1-alpha)*model2
% alpha's can be determined by a subset of tiles
% the best model is found by minimizing the BIC
%exports to R if using a complex, weighted set of rearrangments

%logical: bins x bins x num_param
%e.g. for short events param, denotes are the events in this tile short?
annot_tiles=tiles_annot('length',events,bins,CHR);



% setting up some needed variables
num_param=length(annot_tiles(1,1,:));
alpha=0.5*ones(num_param,1);
if issymmetric(model1)
    model1=triu(model1);
    model1(eye(size(model1))==1)=diag(model1)/2;
end

if issymmetric(model2)
    model2=triu(model2);
    model2(eye(size(model2))==1)=diag(model2)/2;
end

if issymmetric(mfull)
    mfull=triu(mfull);
    mfull(eye(size(mfull))==1)=diag(mfull)/2;
end

%export mfull, model1, model2, annot_tiles to optimize alphas and calculate
%background matrix in R
%Need to do this because we are calculating the log likelihood with the
%continous binomial and the density function for cbinom is in R 
 if complex && weights   
     save('debug_alpha')
     % for 1MB somatic distance cutoff in Xiatong's promiximity matrix
     load('cbinom_alpha.mat');
     %for 50 kb somatic distance cutoff in X's prox matrix 
     %load('2020021050kbcbinom_alpha.mat')   
     disp('loading continous binomial alpha parameters')
 else 
%calculating log factorial--faster than calculating factorial then taking
%the log2
disp('calculating poisson log likelihood approximation')
log_fac(1)=0;
for c1=1:max(mfull(:))
    log_fac(c1+1)=sum(log(1:c1));
end

BIC = zeros(10,1);
alphas = zeros(10,3);

for k1 = 1:10 
nume=sum(mfull(:));
alpha=rand(num_param,1);
nume_1 = sum(mfull(annot_tiles(:,:,1)));
nume_2 = sum(mfull(annot_tiles(:,:,2)));
nume_3 = sum(mfull(annot_tiles(:,:,3)));

% optimize over model
%fincom

%repeat ten times to check for stability
%each row is a new iteration

options = optimoptions('fmincon','DiffMinChange',1e-6,'TolFun',1e-1,'TolX',1e-6);

[opt_alpha, f_bic] = fmincon(@mix_optim_fun,alpha,[],[],[],[],zeros(num_param,1),ones(num_param,1),[],options);

alphas(k1, :) = opt_alpha;
BIC(k1, :) = f_bic;

end 
 
 end 
 
 
 
 
 
 
%set mix_model for the optimal alpha (this is the final mix model that will be returned by the function) 


mix_model=zeros(size(mfull));


for c1=1:num_param
    mix_model(annot_tiles(:,:,c1)) =  mix_model(annot_tiles(:,:,c1)) + opt_alpha(c1)*model1(annot_tiles(:,:,c1))+(1-opt_alpha(c1))*model2(annot_tiles(:,:,c1));
end
%normalize probabilities because multiplying by alpha and (1 - alpha) makes
%the sum of the probabilities slightly less than 1 
mix_model=mix_model/sum(mix_model(:));

    function BIC = mix_optim_fun(alpha)

        % the mix model probability function
        mix_model=zeros(size(mfull));
        for c1=1:num_param
            mix_model(annot_tiles(:,:,c1)) =  mix_model(annot_tiles(:,:,c1)) + alpha(c1)*model1(annot_tiles(:,:,c1))+(1-alpha(c1))*model2(annot_tiles(:,:,c1));
        end
        mix_model=mix_model/sum(mix_model(:));
        
        % the log likelihood function
        % Poisson approximation (faster to compute than binomial and
        % numerically equivalent for large N and small p)
     
        nnz_idc=mix_model>0&mfull>0;
%        nnz_idc=mix_model>0;
         z_idc =  mfull == 0;     

% change loglikelihood to compress long and ic
    mfull_short = mfull(annot_tiles(:, :, 1)); 
    mfull_long = sum(mfull(annot_tiles(:, :, 2)));
    mfull_ic = sum(mfull(annot_tiles(:, :, 3)));
    
    mm_short = mix_model(annot_tiles(:, :, 1));
    mm_long = sum(mix_model(annot_tiles(:, :, 2)));
    mm_ic = sum(mix_model(annot_tiles(:, :, 3)));
    
    %compare mean number of rearrangements with expected value
    m1_short = model1(annot_tiles(:, :, 1));
    m1_long = model1(annot_tiles(:, :, 2));
    m1_ic = model1(annot_tiles(:, :, 3));
    
    m2_short = model2(annot_tiles(:, :, 1));
    m2_long = model2(annot_tiles(:, :, 2));
    m2_ic = model2(annot_tiles(:, :, 3));
    
   % sum(abs(nume*m1_short - mfull_short))
    %sum(abs(nume*m2_short - mfull_short))
    
     sum(nume*m1_long)- mfull_long
    sum(nume*m2_long) - mfull_long
    
    sum(nume*m1_ic) - mfull_ic
    sum(nume*m2_ic) - mfull_ic
     
    
    nnzs_idc = mm_short > 0;
    
    shortLij = sum(mfull_short(nnzs_idc).*log(nume*mm_short(nnzs_idc))-nume*mm_short(nnzs_idc)-log_fac(mfull_short(nnzs_idc)+1)');
    longLij = mfull_long.*log(nume*mm_long)-nume*mm_long-sum(log(1:mfull_long));
   icLij = mfull_ic.*log(nume*mm_ic)-nume*mm_ic-sum(log(1:mfull_ic));
   
   %longLij = binopdf(mfull_long, nume, mm_long);
   %icLij = binopdf(mfull_ic, nume, mm_ic);
   
    %try using binomial distribution instead of poisson approximation
    
         
%mfull are the values of the "poisson" distributed xij, mix_model are the
%pijs, nume is N 
        sLij = sum(sum(mfull(nnz_idc).*log(nume*mix_model(nnz_idc))-nume*mix_model(nnz_idc)-log_fac(mfull(nnz_idc)+1)')) ;
        zLij = sum(sum(-nume*mix_model(z_idc)))  ;    
  
%parameter penality 
 penalty1 = log(nume_1)*(1 - alpha(1))* 8;
 penalty2 = log(nume_2)*(1 - alpha(2)) * 2;
 %penalty3 = (1- alpha(3)) * 0;
 %3 added terms to BIC as penalty
        % the BIC value
    BIC = -2*(sLij + zLij)+log(nume)*num_param
      %BIC = -2*(sLij + zLij)+penalty1 + penalty2 
      
      %BIC with compressed short and ic tiles
      %BIC =  -2 *(shortLij + longLij +icLij) + log(nume)*num_param;
     
      
    end
 
%save('simple_alpha')

 end
%import opt_alpha
        

© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
