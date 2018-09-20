% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function X= gen_newstate_fn(model,Xd,V)

%linear state space equation (CV model)

if ~isnumeric(V)
    if strcmp(V,'noise')
        % V= model.sigma_V*model.B*randn(size(model.B,2),size(Xd,2));
        V= model.G*randn(size(model.G,2),size(Xd,2));
    elseif strcmp(V,'noiseless')
        % V= zeros(size(model.B,1),size(Xd,2));
        V= zeros(size(model.G,1),size(Xd,2));
    end
end

if isempty(Xd)
    X= [];
else
    X= model.F*Xd+ V;
end
