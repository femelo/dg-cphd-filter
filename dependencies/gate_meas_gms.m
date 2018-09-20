% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function z_gate= gate_meas_gms(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2);
plength = size(m,2);

for j=1:plength
    Sj = model.R + model.H*P(:,:,j)*model.H';
    % Changed by Flávio
    % Vs= chol(Sj); det_Sj= prod(diag(Vs))^2; inv_sqrt_Sj= inv(Vs);
    % iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    sqrt_Sj = chol(Sj,'lower');
    nu= z- model.H*repmat(m(:,j),[1 zlength]);
    % dist= sum((inv_sqrt_Sj'*nu).^2);
    dist = sum((sqrt_Sj\nu).^2);
    valid_idx= [valid_idx,find( dist < gamma )];
    
    % Validation based on kd-tree
    % Sj = model.R + model.H*P(:,:,j)*model.H';
    % sqrt_Sj = chol(Sj,'lower');
    % idx_c = rangesearch((sqrt_Sj\z)',(sqrt_Sj\model.H*repmat(m(:,j),[1 zlength]))',sqrt(gamma),'Distance','euclidean','NSMethod','kdtree');
    % idx = idx_c{1,1};
end
valid_idx = unique(valid_idx);
z_gate = z(:,valid_idx);
