function [p, qe, qsolve] = q_solver(R, sij,stat_out)
global firstbin
global lastbin
global intra_chr_a


if abs(sum(R)-2) > 1e-3
    disp('ERROR-solve_generative  R must sum to 2')
end

[nbt, nbt] = size(sij);

qe = R/2; % qe=rand(length(R),1);qe=qe/sum(qe);
maxit=100;
[qsolve,flag1,rr1,iter1] = pcg(sij+eye(size(sij)),R, [],maxit);

if flag1>0
    disp('ERROR in solver convergence');
else
    disp(['Converged on ' num2str(iter1) ' of ' num2str(maxit)]);
end
  
qe = qsolve;
countnegs = find(qsolve<0);  
   
if length(countnegs)>=0 
    disp(['num of negatives in qsolve=' num2str(length(countnegs)) ' min qsolve ' num2str(min(qsolve)) ]);
    qe = projectOntoSimplex(qsolve);
end

p0 = bsxfun(@times,qe,sij);
%p=p0+p0';


sji=sij';

maxit=100;
[qsolve1,flag1,rr1,iter1] = pcg(sji+eye(size(sji)),R, [],maxit);

if flag1>0
    disp('ERROR in solver convergence');
else
    disp(['Converged on ' num2str(iter1) ' of ' num2str(maxit)]);
end
  
qe1 = qsolve1;
countnegs = find(qsolve1<0);  
   
if length(countnegs)>=0 
    disp(['num of negatives in qsolve=' num2str(length(countnegs)) ' min qsolve ' num2str(min(qsolve1)) ]);
    qe1 = projectOntoSimplex(qsolve1);
end
p1 = bsxfun(@times,qe1,sji);


p=p0+p1;

%fix ratio of ic probabilities: non-ic probabilities
%to reflect ratio of ic events:non ic events
%ratio = sum(sum(p(chr_intra, chr_intra)))/sum(sum(p));
ratio = (sum(sum((p(firstbin:lastbin,lastbin+1:end)))) + sum(sum(p(firstbin:lastbin,1:firstbin-1))))/sum(sum(p));


%want ratio to be 1-  intra_chr_a
%p(firstbin:lastbin,lastbin+1:end) = ((1 - intra_chr_a)/ratio)*p(firstbin:lastbin,lastbin+1:end);
%p(firstbin:lastbin,1:firstbin-1) = ((1 - intra_chr_a)/ratio)*p(firstbin:lastbin, 1:firstbin-1);

%normalize
p = 2*p/sum(sum(p));

%recalculate ratio
(sum(sum((p(firstbin:lastbin,lastbin+1:end)))) + sum(sum(p(firstbin:lastbin,1:firstbin-1))))/sum(sum(p))
%check solution
if stat_out

disp('stats on q');
disp('-----------------');
disp(['sum(qsolve) = ' num2str(sum(qsolve))]);
disp(['min qsolve = ' num2str(min(qsolve(:)))]);
disp(['max qsolve = ' num2str(max(qsolve(:)))]);
disp(['mean qsolve = ' num2str(mean(qsolve(:)))]);
disp(['sum(projected-onto-simplex qe) = ' num2str(sum(qe))]);
disp(['min qe = ' num2str(min(qe(:)))]);
disp(['max qe = ' num2str(max(qe(:)))]);
disp(['mean qe = ' num2str(mean(qe(:)))]);
disp(['max diff abs(qe-qsolve) = ' num2str(max(abs(qe(:)-qsolve(:))))]);

disp(['final sol: qsum = ' num2str(sum(qe))]);




%check solution
disp('stats on pij');
disp('-----------------');
disp(['sum(pij) = ' num2str(sum(sum(p0+p0')))]);
disp(['sum(normalized pij) = ' num2str(sum(sum(p)))]);
disp(['min p = ' num2str(min(p(:)))]);
disp(['max p = ' num2str(max(p(:)))]);
disp(['mean p = ' num2str(mean(p(:)))]);
end


