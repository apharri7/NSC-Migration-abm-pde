function [nphi, q, ld, N_KL, tol] = get_KL_coeff(lc, tt, tf, Nmax, cutoff)
% get_KL_coeff constructs a truncated one-dimensional Karhunen-Loeve basis

    global Energy_Cutoff

    if( nargin == 2)
        tf = max(tt); 
    end
    if (nargin <= 3)                 
      [w ld] = eigenvalue_expker(lc, tf, 100);
      Energy_Cutoff = 0.95;
    elseif (nargin == 4 )
      [w ld] = eigenvalue_expker(lc, tf, Nmax);
      Energy_Cutoff = 0.95;    
    elseif (nargin == 5 )
      [w ld] = eigenvalue_expker(lc, tf, Nmax);
      Energy_Cutoff = cutoff;    
      if( cutoff < 0.1 ); Energy_Cutoff = 1 - cutoff;  end 
    else
      disp('wrong number of input'); 
    end  
    
    tol = cumsum(ld)/sum(ld); 
    idx = find(cumsum(ld)/sum(ld)>=Energy_Cutoff);
    N_KL = idx(1);
    ld = ld(1:N_KL);
    
    len=length(tt);
    q = zeros(len, N_KL);
    nphi = zeros(len, N_KL);
    for i=1:N_KL
        q(:,i) = eigenfunction(tt, ld(i), tf, 1/lc);
        nphi(:,i) = q(:,i) * sqrt(ld(i));
    end
end

