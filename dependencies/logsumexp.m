% Written by Daniel S. Bryant

function logsum = logsumexp(w,dim)

%performs log-sum-exp trick to avoid numerical underflow
%input:  w weight vector assumed already log transformed
%output: log(sum(exp(w)))



if nargin == 1
    [val,~] = max(w);
    logsum = log(sum(exp(w-val))) + val;
else
    [val,~] = max(w,[],dim);
    logsum = log(sum(exp(w-val),dim)) + val;
end

