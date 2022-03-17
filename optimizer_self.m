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

function [x, u, new_b, exitflag,info] = optimizer_self(TrackMPC,MPC_vars,ModelParams,n_cars, Y, x, u,  x0, u0,i,convergeflag)

    % reshape X such that it can be used in DP
    X = reshape(x,(MPC_vars.N+1)*ModelParams.nx,1);                       %按照列的顺序进行转换的，也就是第一列读完，读第二列，按列存放（把7行121列的矩阵x重新排列成121*7行一列的矩阵）
    % add other cars to X for DP
    if ~isempty(Y)
        Y = repmat(Y,MPC_vars.N+1,1);                                     %remat将矩阵 A 复制 m×n 块，即把 A 作为 B 的元素，B 由 m×n 个 A 平铺而成  这里是把7行4列的Y重新整理成121行1列的大矩阵
    end
    XX = [X,Y];                                                            %第一列为当前车的信息，第二到五列为障碍车的信息 
    % get boarder points using DP based algorithm
    [new_b,~,grid_isOccupiedBy] = getNewBorders(TrackMPC.traj,TrackMPC.borders,TrackMPC.track_center,XX,n_cars,1,MPC_vars,ModelParams,i,convergeflag);
    
    % formulate MPCC problem and solve it given the DP bundaries
    [xpred, upred,dupred,info] = getMPCmatrices(TrackMPC.traj,MPC_vars,ModelParams,XX,new_b,x,u,x0,u0,grid_isOccupiedBy,i);

%     if (info.exitflag.problem == 0)
%     if (info.exitflag.info == 0) %for yalmip
    if (info.exitflag== 0) 
        exitflag = 0;
        x = xpred;
        u = upred;
    else
        exitflag = 1;
        % x and u stay identical to the initial guess
    end
end