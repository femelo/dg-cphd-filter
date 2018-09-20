% The original code is freely available at http://ba-tuong.vo-au.com/codes.html

function [w_new,x_new,P_new]= gaus_prune(w,x,P,elim_threshold)

idx= find( w > elim_threshold );
w_new= w(idx);
x_new= x(:,idx);
P_new= P(:,:,idx);
