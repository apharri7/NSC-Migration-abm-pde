%
% kernel : exp(- |t_1-t_2|/b) where t \in [0, T] and b : correlation
% length. 
%		N is the length of lambda.
%
%		w is the solution to the equation (w^2-c^2) tan(wT) - 2 c w = 0.
% 	the eigenvalue corresponding to w is 
%			lambda = 2*c / (w^2+c^2)
%

function [w lam] = eigenvalue_expker(b, T, N)

	c = 1/b;      xMax = 100; % if( b < 1 ); xMax = 100;  else; xMax = 1;  end;     
                  dx = xMax/5000; 

	w1 = [0]; 
	maxw1 = 0; maxw2 = 0;
	for x=0:dx:xMax 
			[xzero, fval, exitflag] = fzero(@(x) (x^2-c^2)*tan(x*T) - 2*c*x, x);
			if ( exitflag == 1 && xzero > 0)
					w1 = horzcat(w1, xzero);
					maxw1 = xzero;
			end
	end

	w = unique(sort(w1'));
	w = w(w>eps);
	neww=[]; cnt=1; 
	neww(cnt)=w(1); 
	myeps = 1e-10;
	for n=2:length(w)
		if abs(w(n) - neww(cnt)) > myeps
			neww(cnt+1) = w(n);
			cnt=cnt+1;
		end
	end
			
	w = neww';
	lam = 2*c./(w.^2+c^2); 
    if( length(w)>N )
        w = w(1:N); lam = lam(1:N);
    end

end

