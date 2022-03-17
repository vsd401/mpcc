

function [Xcoefs Ycoefs] = spline5(X,Y,L)
%X Y are coordinates of points (1D), length(X)==length(Y)
%L is the length of each spline. By default L=ones(size(X))
%Xcoefs and Ycoefs are matrices length(X)
assert(length(X)==length(Y),'X and Y should be the same length');


end