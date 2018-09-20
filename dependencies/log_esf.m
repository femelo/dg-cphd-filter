% Written by Flávio Eler De Melo
% Licensed under GNU GPL v3
% Copyright 2018 - Flávio Eler de Melo (flavio.eler@gmail.com)

function log_s = log_esf(log_Z)

%calculate elementary symmetric function using Mahler's recursive formula

if isempty(log_Z)
    log_s= 0;
    return;
end

n_z = length(log_Z);
log_F = zeros(2,n_z);

i_n = 1;
i_nminus = 2;

for n = 1:n_z
    log_F(i_n,1) = logsumexp([log_F(i_nminus,1),log_Z(n)]);
    for k = 2:n
        if k==n
            log_F(i_n,k) = log_Z(n)+log_F(i_nminus,k-1);
        else   
            log_F(i_n,k) = logsumexp([log_F(i_nminus,k),log_Z(n)+log_F(i_nminus,k-1)]);
        end            
    end
    tmp = i_n;
    i_n = i_nminus;
    i_nminus = tmp;
end    

log_s= [0; log_F(i_nminus,:)']; 
