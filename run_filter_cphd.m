function est = run_filter_cphd(model,meas)

% This is the MATLAB code based on the the GMCPHD filter (assuming Poisson clutter) 
% proposed in
% B.-T. Vo, B.-N. Vo and A. Cantoni, "Analytic implementations of the Cardinalized Probability Hypothesis Density Filter," IEEE Trans Signal Processing, Vol. 55,  No. 7, part 2,  pp. 3553-3567, 2007.
% http://ba-ngu.vo-au.com/vo/VVC_CPHD_SP07.pdf
%
% The original code was written by Ba-Tuong Vo and is available at http://ba-tuong.vo-au.com/codes.html
% This code was modified to implement:
% - Pre-computation of factorial factors
% - Re-normalization of intermediate factors to allow for higher number of false alarms
% - Re-normalization of the arguments for computing the elementary symmetric functions
% - Data update using logarithm of factors (to avoid overflow)
%
% Modified by Flávio Eler De Melo
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

filter.N_max= round(2*model.num_targets);      %maximum cardinality number (for cardinality distribution)

filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'silence';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Precompute factors
log_1_n = log(1:filter.N_max);
sum_log_1_n = zeros(1,filter.N_max);
for n = 1:filter.N_max
    sum_log_1_n(n) = sum(log_1_n(1:n));
end

log_lambda_c = log(model.lambda_c);
log_pdf_c = log(model.pdf_c);
log_P_D = log(model.P_D);
log_Q_D = log(model.Q_D);
log_P_S = log(model.P_S);
log_Q_S = log(model.Q_S);

%=== Filtering 

%initial prior
% w_update(1)= eps;
% m_update(:,1)= [0.1;0;0.1;0];
% P_update(:,:,1)= diag([1 1 1 1]).^2;
% L_update = 1;
w_update = [];
m_update = [];
P_update = [];
L_update = 0;
cdn_update= [1; zeros(filter.N_max,1)];

log_min_val = log(eps(0));

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
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            % terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(model.P_S)+(ell-j)*log(model.Q_S))*cdn_update(idxl);
            terms(idxl) = exp(sum_log_1_n(max(ell,1))-sum_log_1_n(max(j,1))-sum_log_1_n(max(ell-j,1))+j*log_P_S+(ell-j)*log_Q_S)*cdn_update(idxl);
        end
        survive_cdn_predict(idxj) = sum(terms);
    end

    %predicted cardinality= convolution of birth and surviving cardinality distribution
    cdn_predict = zeros(filter.N_max+1,1);
    sum_w_birth = sum(model.w_birth);
    for n=0:filter.N_max
        idxn=n+1;
        terms= zeros(filter.N_max+1,1);
        for j=0:n
            idxj= j+1;
            % terms(idxj)= exp(-sum(model.w_birth)+(n-j)*log(sum(model.w_birth))-sum(log(1:n-j)))*survive_cdn_predict(idxj);
            terms(idxj)= exp(-sum_w_birth+(n-j)*log(sum_w_birth)-sum_log_1_n(max(1,n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
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
   
    t_start = tic; 
    %pre calculation for Kalman update parameters
    if m~=0
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
    end
    
    %pre calculation for elementary symmetric functions
    XI_vals = zeros(m,1);                        %arguments to esf
    for ell=1:m
       XI_vals(ell) = model.P_D*w_predict'*qz_temp(:,ell)/model.pdf_c;
    end
    
    sum_w_predict = sum(w_predict);
    den_factor = sum_w_predict * model.lambda_c;
    esfvals_E = esf(XI_vals / den_factor);      %calculate esf for entire observation set
    esfvals_D = zeros(m,m);                     %calculate esf with each observation index removed one-by-one
    for ell=1:m
        esfvals_D(:,ell) = esf([XI_vals(1:ell-1);XI_vals(ell+1:m)] / den_factor);
    end
    
    %pre calculation for upsilons
    % upsilon0_E = zeros(filter.N_max+1,1);
    % upsilon1_E = zeros(filter.N_max+1,1);
    % upsilon1_D = zeros(filter.N_max+1,m);
    log_upsilon0_E = log_min_val*ones(filter.N_max+1,1);
    log_upsilon1_E = log_min_val*ones(filter.N_max+1,1);
    log_upsilon1_D = log_min_val*ones(filter.N_max+1,m);
    
    log_sum_w_predict = log(sum(w_predict));
    log_esfvals_E = log(esfvals_E);
    log_esfvals_D = log(esfvals_D);
    for n=0:filter.N_max
        idxn= n+1;
        
        % terms0_E = zeros(min(m,n)+1,1);  %calculate upsilon0_E(idxn)
        log_terms0_E= log_min_val*ones(min(m,n)+1,1);
        for j=0:min(m,n)
            idxj= j+1;
            % Tweaked by Flávio
            % terms0_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-j))+(n-j)*log(model.Q_D)-j*log(sum(w_predict)))*esfvals_E(idxj);
            % terms0_E(idxj) = exp(-model.lambda_c+(-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-j))+(n-j)*log(model.Q_D)-j*log(sum(w_predict)))*esfvals_E(idxj);
            % terms0_E(idxj) = exp(-model.lambda_c+(-j)*log_lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-j,1))+(n-j)*log_Q_D-j*log_sum_w_predict)*esfvals_E(idxj);
            log_terms0_E(idxj) = -model.lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-j,1))+(n-j)*log_Q_D+log_esfvals_E(idxj);
        end
        % upsilon0_E(idxn)= sum(terms0_E);
        log_upsilon0_E(idxn)= logsumexp(log_terms0_E);
        
        % terms1_E= zeros(min(m,n)+1,1);  %calculate upsilon1_E(idxn)
        log_terms1_E= log_min_val*ones(min(m,n)+1,1);
        for j=0:min(m,n)
            idxj= j+1;
            if n>=j+1
                % terms1_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_E(idxj);
                % terms1_E(idxj) = exp(-model.lambda_c+(-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_E(idxj);
                % terms1_E(idxj) = exp(-model.lambda_c+(-j)*log_lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-(j+1),1))+(n-(j+1))*log_Q_D-(j+1)*log_sum_w_predict)*esfvals_E(idxj);
                log_terms1_E(idxj) = -model.lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-(j+1),1))+(n-(j+1))*log_Q_D-log_sum_w_predict+log_esfvals_E(idxj);
            end
        end
        % upsilon1_E(idxn)= sum(terms1_E);
        log_upsilon1_E(idxn)= logsumexp(log_terms1_E);
        
        if m~= 0 %calculate upsilon1_D(idxn,:) if m>0
            % terms1_D= zeros(min((m-1),n)+1,m);
            log_terms1_D= log_min_val*ones(min((m-1),n)+1,m);
            for ell=1:m
                for j=0:min((m-1),n)
                    idxj= j+1;
                    if n>=j+1
                        % terms1_D(idxj,ell) = exp(-model.lambda_c+((m-1)-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_D(idxj,ell);
                        % terms1_D(idxj,ell) = exp(-model.lambda_c+((-1)-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_D(idxj,ell);
                        % terms1_D(idxj,ell) = exp(-model.lambda_c+((-1)-j)*log_lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-(j+1),1))+(n-(j+1))*log_Q_D-(j+1)*log_sum_w_predict)*esfvals_D(idxj,ell);
                        log_terms1_D(idxj,ell) = -model.lambda_c-log_lambda_c+sum_log_1_n(max(n,1))-sum_log_1_n(max(n-(j+1),1))+(n-(j+1))*log_Q_D-log_sum_w_predict+log_esfvals_D(idxj,ell);
                    end
                end
            end
            % upsilon1_D(idxn,:)= sum(terms1_D,1);
            log_upsilon1_D(idxn,:)= logsumexp(log_terms1_D,1);
        end
    end
    

    %missed detection term 
    % w_update = (upsilon1_E'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.Q_D*w_predict;
    log_cdn_predict = log(cdn_predict);
    log_w_predict = log(w_predict);
    log_w_update = logsumexp(log_upsilon1_E+log_cdn_predict) - logsumexp(log_upsilon0_E+log_cdn_predict) + log_Q_D + log_w_predict;
    m_update = m_predict;
    P_update = P_predict;
    
    if m~=0
        %m detection terms 
        %Kalman update precalculated and stored
        for ell=1:m
            % w_temp = (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.P_D.*qz_temp(:,ell)/model.pdf_c.*w_predict;
            % w_update = cat(1,w_update,w_temp);
            log_w_temp = logsumexp(log_upsilon1_D(:,ell)+log_cdn_predict) -logsumexp(log_upsilon0_E+log_cdn_predict) +log_P_D +log(qz_temp(:,ell)) -log_pdf_c +log_w_predict;
            log_w_update = cat(1,log_w_update,log_w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end    
    w_update = exp(log_w_update);
    
    %---cardinality update
    % cdn_update= upsilon0_E.*cdn_predict;
    log_cdn_update = log_upsilon0_E+log_cdn_predict;
    cdn_update= exp(log_cdn_update-max(log_cdn_update));
    cdn_update= cdn_update/sum(cdn_update);
            
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
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    est.N(k) = min(length(w_update),map_cdn);
    est.S(k) = ([0:filter.N_max].^2)*cdn_update-([0:filter.N_max]*cdn_update)^2;
    [~,idx_m_srt]= sort(-w_update);
    est.X{k} = m_update(:,idx_m_srt(1:est.N(k)));
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #tot phd=' num2str(sum(w_update),4),...
         ' #eap cdn=' num2str([0:filter.N_max]*cdn_update,4),...
         ' #var cdn=' num2str([0:filter.N_max].^2*cdn_update-([0:filter.N_max]*cdn_update)^2,4),...
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
            
