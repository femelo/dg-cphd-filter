% The original code is freely available at http://ba-tuong.vo-au.com/codes.html
% Modified by Flávio Eler De Melo (flavio.eler@gmail.com)
function model = gen_model(P_D,lambda_c,N)

model.num_time_steps = 100;

% basic parameters
model.x_dim= 4;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector

% dynamical model parameters (CV model)
model.T= 1;                                           %sampling period
model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
model.F= [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
model.sigma_v = 5;
model.Q0 = [model.T^3/3 model.T^2/2; model.T^2/2 model.T];
model.Q = (model.sigma_v^2)*blkdiag(model.Q0,model.Q0);
model.G = chol(model.Q,'lower');

% survival/death parameters
model.P_S = 0.98;
model.Q_S= 1-model.P_S;

% birth parameters (Poisson birth model, multiple Gaussian components)
model.L_birth= 5;                                                     %no. of Gaussian birth terms
model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
model.m_birth= zeros(model.x_dim,model.L_birth);                      %means of Gaussian birth terms 
model.B_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %std of Gaussian birth terms
model.P_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %cov of Gaussian birth terms

lambda_B = (N+5)/200;

model.w_birth(1)= lambda_B/5;                                                 %birth term 1
model.m_birth(:,1)= [ -500; 0; -500; 0 ];
model.B_birth(:,:,1)= diag([ 500/2.5; 10/2.5; 500/2.5; 10/2.5 ]);
model.P_birth(:,:,1)= model.B_birth(:,:,1)*model.B_birth(:,:,1)';
    
model.w_birth(2)= lambda_B/5;                                                 %birth term 2
model.m_birth(:,2)= [ -500; 0; +500; 0 ];
model.B_birth(:,:,2)= diag([ 500/2.5; 10/2.5; 500/2.5; 10/2.5 ]);
model.P_birth(:,:,2)= model.B_birth(:,:,2)*model.B_birth(:,:,2)';

model.w_birth(3)= lambda_B/5;                                                 %birth term 3
model.m_birth(:,3)= [ +500; 0; -500; 0 ];
model.B_birth(:,:,3)= diag([ 500/2.5; 10/2.5; 500/2.5; 10/2.5 ]);
model.P_birth(:,:,3)= model.B_birth(:,:,3)*model.B_birth(:,:,3)';

model.w_birth(4)= lambda_B/5;                                                 %birth term 4
model.m_birth(:,4)= [ +500; 0; +500; 0 ];
model.B_birth(:,:,4)= diag([ 500/2.5; 10/2.5; 500/2.5; 10/2.5 ]);
model.P_birth(:,:,4)= model.B_birth(:,:,4)*model.B_birth(:,:,4)';

model.w_birth(5)= lambda_B/5;                                                 %birth term 5
model.m_birth(:,5)= [ +0; 0; +0; 0 ];
model.B_birth(:,:,5)= diag([ 500/2.5; 10/2.5; 500/2.5; 10/2.5 ]);
model.P_birth(:,:,5)= model.B_birth(:,:,5)*model.B_birth(:,:,5)';

% observation model parameters (noisy x/y only)
model.H= [ 1 0 0 0 ; 0 0 1 0 ];    %observation matrix
model.D= diag([ 5; 5 ]); 
model.R= model.D*model.D';         %observation noise covariance

% detection parameters
model.P_D= P_D;         %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.lambda_c= lambda_c;                                   %poisson average rate of uniform clutter (per scan)
model.range_c= [ -1000 1000; -1000 1000 ];                  %uniform clutter region
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
model.num_targets = N;
