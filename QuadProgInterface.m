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

function [ X,U,dU,info ] = QuadProgInterface(stage,MPC_vars,ModelParams)
nx = ModelParams.nx;
nu = ModelParams.nu;
ng = 2;
ncarid=4;
nz = nx+2*nu;
nxu = nx+nu;
N = MPC_vars.N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = zeros(nz*N+nxu,nz*N+nxu);
f = zeros(nz*N+nxu,1);
for i = 1:N+1
    if i<N+1
        H_i = blkdiag(stage(i).Qk,stage(i).Rk);                           %所以最终H的组成为【轨迹误差ec,el,    三个控制量的权重    三个控制量增量的权重】12/X*H*X，这里的X是1*13维的量，分别是【[X Y PHI Vx Vy omega 轨迹SITA] 三个控制量[电机控制信号  方向盘转角 车速] 三个控制增量[电机控制信号  方向盘转角 车速]】
        H((i-1)*nz+[1:nz],(i-1)*nz+[1:nz]) = H_i;
    
        f((i-1)*nz+[1:nxu]) = stage(i).fk;
    else
        H((i-1)*nz+[1:nxu],(i-1)*nz+[1:nxu]) = stage(i).Qk;
    
        f((i-1)*nz+[1:nxu]) = stage(i).fk;
    end
end

H = 0.5*(H+H');

Aeq = zeros(nxu*(N+1),nz*N+nxu);
beq = zeros(nxu*(N+1),1);

x0scale = blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0];

Aeq(1:nxu,1:nxu) = eye(nxu);
beq(1:nxu) = x0scale;


for i = 1:N
%     if i == 1
%         Aeq((i-1)*nxu+[1:nxu],(i-1)*nz+[1:2*nz]) = [zeros(nxu,nxu),stage(i).Bk,-eye(nxu),zeros(nxu,nu)];
%         beq((i-1)*nxu+[1:nxu]) = -stage(i).Ak*x0scale-stage(i).gk;
    if i<N
        Aeq((i)*nxu+[1:nxu],(i-1)*nz+[1:2*nz]) = [-stage(i).Ak,-stage(i).Bk,eye(nxu),zeros(nxu,nu)];
        beq((i)*nxu+[1:nxu]) = stage(i).gk;
    else
        Aeq((i)*nxu+[1:nxu],(i-1)*nz+[1:(nz+nxu)]) = [-stage(i).Ak,-stage(i).Bk,eye(nxu)];
        beq((i)*nxu+[1:nxu]) = stage(i).gk;
    end
end

% A = [];%ones(1,nz*N+nxu);
% b = [];%zeros(ng*N,1);

A = zeros((ng+ncarid+1)*N,nz*N+nxu);
b = zeros((ng+ncarid+1)*N,1);

for i = 1:N
    A((i-1)*ng+1,i*nz+[1:nxu]) = stage(i+1).Ck;
    b((i-1)*ng+1) = stage(i+1).ug;
    
    A((i-1)*ng+2,i*nz+[1:nxu]) = -stage(i+1).Ck;
    b((i-1)*ng+2) = -stage(i+1).lg;
end 
for i=1:N
    A(ng*N+(i-1)*ncarid+[1:ncarid],i*nz+[1:nxu]) = stage(i+1).Dk;  
    b(ng*N+(i-1)*ncarid+[1:ncarid]) = stage(i+1).dk;
end 

for i=1:N
    A((ng+ncarid)*N,i*nz+[1:nxu]) = stage(i+1).Ek;  
    b((ng+ncarid)*N+i) = stage(i+1).ek;
end 

LB = zeros(nz*N+nxu,1);
UB = zeros(nz*N+nxu,1);

for i=1:N+1
    if i<N+1
        LB((i-1)*nz+[1:nz]) = stage(i).lb;
        UB((i-1)*nz+[1:nz]) = stage(i).ub;
    else
        LB((i-1)*nz+[1:nxu]) = stage(i).lb(1:nxu);
        UB((i-1)*nz+[1:nxu]) = stage(i).ub(1:nxu);
    end
end
% options = optimoptions(@quadprog,'MaxIterations',100,'Display','off');
% options = optimoptions('quadprog','Algorithm','active-set','Display','off');
tic
[z,~,exitflag] = quadprog(H,f,A,b,Aeq,beq,LB,UB,[]);%,options  
% [z,~,exitflag] = quadprog(H,f,A,b,Aeq,beq,LB,UB,[],options);
QPtime = toc;

X = zeros(nx,N+1);
U = zeros(nu,N); 
dU = zeros(nu,N);

if exitflag == 1 
    for i = 1:N+1
        X(1:nx,i) = MPC_vars.invTx*z((i-1)*nz+[1:nx]);                    %取z（1:7）为X状态量

        if i>1
            U(1:nu,i-1) = MPC_vars.invTu*z((i-1)*nz+nx+[1:nu]);           %取z（7:10）为U
        end

        if i<=N
            dU(1:nu,i) = z((i-1)*nz+nxu+[1:nu]);                          %取z（11:13）为det_U
        end    
    end
end

info.QPtime = QPtime;
if exitflag == 1
    info.exitflag = 0;
else
    info.exitflag = 1;
end

end
