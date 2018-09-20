% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function [qz_update,m_update,P_update] = kalman_update_multiple(z,model,m,P,log_flag)

if nargin == 4
    log_flag = false;
end
plength= size(m,2);
zlength= size(z,2);

qz_update= zeros(plength,zlength);
m_update = zeros(model.x_dim,plength,zlength);
P_update = zeros(model.x_dim,model.x_dim,plength);

if log_flag
    for idxp=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single_log(z,model.H,model.R,m(:,idxp),P(:,:,idxp));
        qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
    end
else
    for idxp=1:plength
        [qz_temp,m_temp,P_temp] = kalman_update_single(z,model.H,model.R,m(:,idxp),P(:,:,idxp));
        qz_update(idxp,:)   = qz_temp;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
    end
end

function [qz_temp,m_temp,P_temp] = kalman_update_single(z,H,R,m,P)

mu = H*m;
S  = R+H*P*H';
Vs = chol(S);
det_S= prod(diag(Vs))^2;
% inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K  = P*H'/S;

% qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))))';
nu = z-repmat(mu,[1 size(z,2)]);
qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*sum((Vs\nu).^2,1))';
m_temp = repmat(m,[1 size(z,2)]) + K*nu;
P_temp = P-K*H*P;

function [qz_temp,m_temp,P_temp] = kalman_update_single_log(z,H,R,m,P)

mu = H*m;
S  = R+H*P*H'; Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K  = P*H'*iS;

qz_temp = (-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(mu,[1 size(z,2)]),iS*(z-repmat(mu,[1 size(z,2)]))))';
m_temp = repmat(m,[1 size(z,2)]) + K*(z-repmat(mu,[1 size(z,2)]));
P_temp = (eye(size(P))-K*H)*P;


