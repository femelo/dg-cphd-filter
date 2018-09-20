function est = run_filter_dgcphd(model,meas)

% This is the MATLAB code for the Gaussian Mixture implementation of the 
% Discrete-Gamma CPHD filter proposed in
% F. E. De Melo and S. Maskell, "A CPHD approximation based on a discrete-Gamma cardinality model," IEEE Trans Signal Processing
%
% Licensed under GNU GPL v3
% Copyright 2018 - Flávio Eler de Melo (flavio.eler@gmail.com)
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.S= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.L_max= 300;                  %limit on number of Gaussians
filter.elim_threshold= 1e-5;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold

filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'silence';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
L_update = 0;
w_update = [];
m_update = [];
P_update = [];

a_update = 0;
b_update = 1;

N_update = a_update/b_update;

P_S = model.P_S;
Q_S = model.Q_S;
P_D = model.P_D;
Q_D = model.Q_D;
log_P_D = log(model.P_D);
log_Q_D = log(model.Q_D);
log_P_S = log(model.P_S);
log_Q_S = log(model.Q_S);
pdf_c = model.pdf_c;
log_pdf_c = log(pdf_c);
lambda_c = model.lambda_c;
log_lambda_c = log(lambda_c);

prd_time = 0;
gat_time = 0;
upd_time = 0;
mgm_time = 0;
%recursive filtering
for k=1:meas.K
    t_start = tic;    
    %---intensity prediction 
    [m_predict,P_predict] = kalman_predict_multiple(model,m_update,P_update);                       %surviving components
    w_predict= model.P_S*w_update;                                                                  %surviving weights

    m_predict= cat(2,model.m_birth,m_predict); P_predict=cat(3,model.P_birth,P_predict);            %append birth components
    w_predict= cat(1,model.w_birth,w_predict);                                                      %append birth weights
                                                 
    L_predict= model.L_birth+L_update;                                                              %number of predicted components
    
    % Cardinality prediction
    % Predicted number of targets
    N_predict = sum(model.w_birth) +model.P_S*a_update/b_update;
    % Predicted variance on number of targets
    Sn_predict = N_predict +(model.P_S^2)*N_update*(1/b_update -1);
    
    % Predicted parameters of the binomial distribution
    a_predict = (N_predict^2)/Sn_predict;
    b_predict = N_predict/Sn_predict;
    prd_time = prd_time +toc(t_start);

    %---gating
    t_start = tic;
    if filter.gate_flag
        meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict,P_predict);        
    end
    gat_time = gat_time +toc(t_start);
  
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %pre calculation for Kalman update parameters
    t_start = tic;
    if m~=0
        % [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
        [log_qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict,true);
    end
    
    %pre calculation for elementary symmetric functions
    % XI_vals = zeros(1,m);
    % for j=1:m
    %    XI_vals(j) = model.P_D*w_predict'*qz_temp(:,j)/model.pdf_c;
    % end
    log_w_predict = log(w_predict);
    log_factor_pred = log(a_predict)-log(b_predict);
    factor_pred = a_predict/b_predict;
    XI_vals = zeros(1,m);
    for j=1:m
       XI_vals(j) = exp( log_P_D +logsumexp( log_qz_temp(:,j)+log_w_predict ) -log_pdf_c) + eps(0);
    end
    
    % Evaluate polynomials
    % Theta = poly_theta(m+2,model.Q_D,a_predict,b_predict,N_predict);
    log_Theta = poly_theta(m+2,model.Q_D,a_predict,b_predict,N_predict);
    
    % Evaluate 
    % v = (XI_vals/((a_predict/b_predict)*lambda_c));
    % e_s = esf(v)';
    v = (XI_vals/(factor_pred*lambda_c));
    log_e_s = log(esf(v)');
    
    offset = 1;
    % terms_num_q = Theta(offset+(1:m+1)).*e_s;
    % terms_den_q = Theta(offset+(0:m  )).*e_s;
    log_terms_num_q = log_Theta(offset+(1:m+1))+log_e_s;
    log_terms_den_q = log_Theta(offset+(0:m  ))+log_e_s;
    
    % L_q = sum(terms_num_q)/sum(terms_den_q);
    log_L_q = logsumexp(log_terms_num_q)-logsumexp(log_terms_den_q);
    
    % Coefficients for computing cardinality moments
    % terms_num_p = (0:m).*terms_den_q;
    % l_p =  sum(terms_num_p)/sum(terms_den_q);
    
    % terms_num_r = ((0:m).^2).*terms_den_q;
    % l_r =  sum(terms_num_r)/sum(terms_den_q);
    
    % terms_num_q2 = (0:m).*terms_num_q;
    % l_q = sum(terms_num_q2)/sum(terms_den_q);
    
    % terms_num_s = Theta(offset+(2:m+2)).*e_s;
    % l_s = sum(terms_num_s)/sum(terms_den_q);

    % Compute coefficients using log-factors
    log_0_m = log((0:m)+eps(0));
    log_terms_num_p = log_0_m+log_terms_den_q;
    log_l_p =  logsumexp(log_terms_num_p)-logsumexp(log_terms_den_q);
    
    log_terms_num_r = 2*log_0_m+log_terms_den_q;
    log_l_r =  logsumexp(log_terms_num_r)-logsumexp(log_terms_den_q);
    
    log_terms_num_q2 = log_0_m+log_terms_num_q;
    log_l_q = logsumexp(log_terms_num_q2)-logsumexp(log_terms_den_q);
    
    log_terms_num_s = log_Theta(offset+(2:m+2))+log_e_s;
    log_l_s = logsumexp(log_terms_num_s)-logsumexp(log_terms_den_q);
        
    % L_p = zeros(m,1);
    log_L_p = zeros(m,1)+log(eps(0));
    for ell=1:m
        % v = [XI_vals(1:ell-1), XI_vals(ell+1:m)]/...
        %     ((a_predict/b_predict)*lambda_c);
        
        % e_s = esf(v)';
        
        % terms_num = Theta(offset+(1:m)).*e_s;
        % L_p(ell) = sum(terms_num)/sum(terms_den_q);
        
        % With log-factors
        v = [XI_vals(1:ell-1), XI_vals(ell+1:m)]/...
            (factor_pred*lambda_c);
        
        log_e_s = log(esf(v))';
        
        log_terms_num = log_Theta(offset+(1:m))+log_e_s;
        log_L_p(ell) = logsumexp(log_terms_num)-logsumexp(log_terms_den_q);
    end
        
    %missed detection term 
    % w_update = (L_q/(a_predict/b_predict))*model.Q_D*w_predict;
    log_w_update = log_L_q-log_factor_pred+log_Q_D+log_w_predict;
    m_update = m_predict;
    P_update = P_predict;
    
    % Normalizing term
    % norm_term = lambda_c;
    
    if m~=0
        %m detection terms 
        %Kalman update precalculated and stored
        for ell=1:m
            % w_temp = (L_p(ell)/(a_predict/b_predict))*(model.P_D*(qz_temp(:,ell)/model.pdf_c)/norm_term).*w_predict;
            % w_update = cat(1,w_update,w_temp);
            log_w_temp = log_L_p(ell)-log_factor_pred+log_P_D+log_qz_temp(:,ell)-log_pdf_c-log_lambda_c+log_w_predict;
            log_w_update = cat(1,log_w_update,log_w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end
    w_update = exp(log_w_update);
    
    % Cardinality update
    % Updated number of targets
    % N_update = ( l_p +L_q*model.Q_D );
    N_update = exp(log_l_p) +exp(log_L_q+log_Q_D);
    % Updated variance on number of targets
    % Sn_update = (l_r-l_p) +2*l_q*model.Q_D +l_s*model.Q_D^2 -N_update^2 +N_update;
    Sn_update = exp(log_l_r)-exp(log_l_p) +2*exp(log_l_q+log_Q_D) +exp(log_l_s+2*log_Q_D) -N_update^2 +N_update;
    
    % Updated parameters of the discrete gamma distribution    
    % a_update = round( (N_update^2)/Sn_update );
    a_update = (N_update^2)/Sn_update;
    b_update = N_update/Sn_update;
    
    %---mixture management
    L_posterior= length(w_update);
    upd_time = upd_time +toc(t_start);
    
    %pruning, merging, capping
    t_start = tic;
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);               L_cap  = length(w_update);
    mgm_time = mgm_time +toc(t_start);

    L_update= L_cap;
    
    %--- state extraction
    % est.N(k) = min(length(w_update),map_cdn);
    est.N(k) = round(N_update);
    est.S(k) = Sn_update;
    [~,idx_m_srt]= sort(-w_update);
    est.X{k} = m_update(:,idx_m_srt(1:round(est.N(k))));
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #tot phd=' num2str(sum(w_update),4),...
         ' #eap cdn=' num2str(N_update,4),...
         ' #var cdn=' num2str(Sn_update,4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);
    end
end
est.prd_time = prd_time;
est.gat_time = gat_time;
est.upd_time = upd_time;
est.mgm_time = mgm_time;

end

% function [Theta] = poly_theta(r,s,alpha,beta,N)
function [log_Theta] = poly_theta(r,s,alpha,beta,N)
    log_s = log(s);
    log_z = -beta+log_s;
    % Psi = zeros(1,r+1);
    log_Psi = zeros(1,r+1)+log(eps(0));
    
    % Number of terms to approximate    
    epsilon = realmin;
    nu = N;
    for k = 1:30
        nu = nu -((alpha-1)*log(nu)-beta*nu-log(epsilon))/(2*((alpha-1)/nu -beta));
    end
    lb = 1;
    % ub = round(real(nu));
    ub = max(round(real(nu)),r);
    
    n = lb:ub;
    n_log_z = n*log_z;
    log_n = log(n);
    am1_log_n = (alpha-1)*log_n;
    theta = am1_log_n +n_log_z;
    theta_max = max(theta);
    theta_min = min(theta);
    % d_theta = theta-theta_max;
    d_theta = theta-(theta_min+theta_max)/2;
    
    % Psi(1) = sum(exp(d_theta));
    log_Psi(1) = logsumexp(d_theta);
    
    l = 1:r;
    
    % Different method
    log_s_j = log_s;
    log_n_j = log_n;
    for j = l
        % Psi(j+1) = sum(real(exp(d_theta+log_n_j-log_s_j)));
        log_Psi(j+1) = logsumexp(d_theta+log_n_j-log_s_j);
        % log_n_j = log_n_j+log(n-j);
        log_n_j = log_n_j+log(max(n-j,0));
        log_s_j = log_s_j+log_s;
    end
    
    % Theta = Psi/max(abs(Psi));
    [~,i_max] = max(abs(log_Psi));
    log_Theta = log_Psi-log_Psi(i_max);
end


            
