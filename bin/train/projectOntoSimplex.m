function x=projectOntoSimplex(y)
  
  % from http://arxiv.org/pdf/1309.1541v1.pdf "projection onto a probability simplex"
  % following their notation:
  D = length(y);
  x=0*y;

  [u,si] = sort(y,'descend');
  %compute rho
  rho=1;
  for j=1:D
    rarg(j) = u(j)+1/j*(1-sum(u(1:j)));
    if rarg(j)>0
       rho=j;
    end
  end
  lambda=1/rho*(1-sum(u(1:rho)));
  for j=1:D
    x(j)=max(y(j)+lambda,0);
  end
  
  %end of method 3: true projection onto simplex method.

  
