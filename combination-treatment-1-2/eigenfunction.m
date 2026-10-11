%
% The eigenfunction corresponding to eigenvalue lambda.
% INPUT:
%   tf = final time
%   c = 1/cor
%

function f = eigenfunction(t,lambda,tf,c)
    f = phibar(t, sqrt(2*c/lambda-c^2), tf, c);
end
