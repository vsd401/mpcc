% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dppx]=getSplineDerivatives(ppx)
%[dppx]=getSplineDerivatives(ppx)
%Compute the descriptor of the derivative of a spline
%ppx descriptor of the spline
if(ppx.order<1)
    error('order of spline should be more than 1');
end
dppx=ppx;

dppx.coefs=ppx.coefs(:,1:end-1).*kron(ones(size(ppx.coefs,1),1),(size(ppx.coefs,2)-1):-1:1); %%对多项式求导
dppx.order = ppx.order-1;
end
