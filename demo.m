% This is a demo script for the DG-CPHD filter proposed in
% F. E. De Melo and S. Maskell, "A CPHD approximation based on a discrete-Gamma cardinality model," IEEE Trans Signal Processing
%
% Licensed under GNU GPL v3
% Copyright 2018 - Flávio Eler de Melo (flavio.eler@gmail.com)
%---

addpath('./dependencies/')
addpath('./dependencies/export_fig/')
rng(1);
p_d = 0.95;  % probability of detection
lambda = 50; % expected number of false alarms per frame
N_t = 30;    % number of targets
model= gen_model(p_d,lambda,N_t);
truth= gen_truth(model);
meas=  gen_meas(model,truth);

% Run sequence
% 1: PHD
% 2: CPHD
% 3: DG-PHD
seq = [1 2 3];
lbl = {'phd', 'cphd', 'dgcphd'};
txt = {'PHD filter', 'CPHD filter', 'DG-CPHD filter'};

% profile on;
% for j = seq
%     t_0 = tic;
%     result.(lbl{j}).est = feval(['run_filter_' lbl{j}],model,meas);
%     result.(lbl{j}).time = toc(t_0)-result.(lbl{j}).est.m_time;
% end
%% profile viewer;

%profile on;
for j = seq
    t_0 = tic;
    result.(lbl{j}).est = feval(['run_filter_' lbl{j}],model,meas);
    result.(lbl{j}).time = toc(t_0)-result.(lbl{j}).est.mgm_time;
end
%profile viewer;
%%
fprintf('\n   PHD: %05.2f seconds\n',result.phd.time);
fprintf('\n  CPHD: %05.2f seconds\n',result.cphd.time);
fprintf('\nDG-PHD: %05.2f seconds\n',result.dgcphd.time);

% fprintf('\n    PHD: %05.2f seconds\n   CPHD: %05.2f seconds\nDG-CPHD: %05.2f seconds\n',result.phd.time,result.cphd.time,result.dgcphd.time);

generate_plots(model,truth,meas,result.phd.est,result.cphd.est,result.dgcphd.est);

