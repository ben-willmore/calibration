function y=jdecimate(x);
% function y=jdecimate(x);
%   jan's decimation function. 
%   Returns a y of half the length of x, where y(ii,:)= x(2*ii-1,:)+x(2*ii),:) / 2
%   i.e. it returns the row-wise average between successive pairs of values in rows of x.

xlen=size(x,1);
evenIdx=2:2:xlen;
oddIdx=evenIdx-1;
y=x(oddIdx,:)+x(evenIdx,:);
y=y/2;