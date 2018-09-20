% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function Xc= get_comps(X,c)

if isempty(X)
    Xc= [];
else
    Xc= X(c,:);
end
