
% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the 'License');
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an 'AS IS' BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,U,dU,info] = getMPCmatrices(traj,MPC_vars,ModelParams,X_all,borders,Xhor,Uhor,x0,u0,grid_isOccupiedBy,SimN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each stage in the horizon compute the necessary system and cost %%%%%
% matricies and safe them in a array of structs 'stage' %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost scaling for numerics
costScale = MPC_vars.costScale;
% init stage struct
stage = struct([]);
%% Generate MPCC problem
% initial state (including previus input for input rate cost)

% state is augmented by the inputs and rate inputs are introduced to
% formulate rate input costs and constraints while retaining a block
% spares formulation
% given x_k+1 = A x_k + B u_k
% do the following state augmentation
% s_k = [x_k,u_k-1], v_k = du_k
% with the following linear system
% s_k+1 = [A B;0 I] s_k + [B;I] v_k

stage(1).x0 = x0;                                                         %状态量
stage(1).u0 = u0;                                                         %控制量

for i = 1:MPC_vars.N
    Xk = Xhor(:,i);
    Uk = Uhor(:,i);
    % generate quadratic state(-input) cost
    stage(i).Qk = costScale*generateH(traj,MPC_vars,ModelParams,Xk,i);
    % quadratic rate input cost 
    stage(i).Rk = costScale*2*diag([MPC_vars.rdD; MPC_vars.rdDelta; MPC_vars.rdVtheta]);       %控制量增量的权重Rk    Qk和Rk都是为了计算目标函数H的 
    % linear state(-input) cost
    stage(i).fk = costScale*generatef(traj,MPC_vars,ModelParams,Xk,i);                         %
    % linearized dynamics
    [stage(i).Ak,stage(i).Bk,stage(i).gk] = getEqualityConstraints(Xk,Uk,MPC_vars,ModelParams);%等式约束，也就是动力学方程
    % linearized track constraints
    [stage(i).Ck, stage(i).ug, stage(i).lg] = getInequalityConstraints(borders(max(i-1,1),:),MPC_vars,ModelParams);%不等式约束，也就是道路边界约束，三角形定理
    [stage(i).Dk, stage(i).dk] = getInequalityConstraints_OA(MPC_vars,ModelParams,grid_isOccupiedBy,X_all((i-1)*ModelParams.nx+1:i*ModelParams.nx,:),SimN);
    [stage(i).Ek, stage(i).ek] = getInequalityConstraints_OB(traj,MPC_vars,ModelParams,Xk);
    % bounds
    [stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);                               %对优化变量的约束   即：-10000<X<10000 【1, 1,10】为U的约束；[0.25 , 0.1, 10]为det_U的约束
end
% terminal stage
i = MPC_vars.N+1; 
Xk = Xhor(:,i);
% generate quadratic state(-input) cost
stage(i).Qk = costScale*generateH(traj,MPC_vars,ModelParams,Xk,i);
% quadratic rate input cost 
stage(i).Rk = costScale*2*diag([MPC_vars.rdD; MPC_vars.rdDelta; MPC_vars.rdVtheta]);
% linear state(-input) cost
stage(i).fk = costScale*generatef(traj,MPC_vars,ModelParams,Xk,i);
% linearized track constraints
[stage(i).Ck, stage(i).ug, stage(i).lg] = getInequalityConstraints(borders(i-1,:),MPC_vars,ModelParams);
[stage(i).Dk, stage(i).dk] =  getInequalityConstraints_OA(MPC_vars,ModelParams,grid_isOccupiedBy,X_all((i-1)*ModelParams.nx+1:i*ModelParams.nx,:),SimN);
[stage(i).Ek, stage(i).ek] = getInequalityConstraints_OB(traj,MPC_vars,ModelParams,Xk);
% bounds
[stage(i).lb, stage(i).ub] = getBounds(MPC_vars,ModelParams);
%% Call solver interface
if strcmp(MPC_vars.interface,'Yalmip') 
    % yalmip based interface (very slow)
    [X,U,dU,info] = YalmipInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface,'CVX')
    % CVX based interface (slow)
    [X,U,dU,info] = CVXInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface,'hpipm')
    % hpipm interface (prefered)
    [X,U,dU,info] = hpipmInterface(stage,MPC_vars,ModelParams);
elseif strcmp(MPC_vars.interface,'quadprog')
    % quadprog interface (replace quadprog with a better solver if possible)
    [X,U,dU,info] = QuadProgInterface(stage,MPC_vars,ModelParams);
else
    error('invalid optimization interface')
end

end

% GENERATING Q
function Qk = generateH(pathinfo,MPC_vars,ModelParams,Xk,i)
    % get linearized contouring and lag errors
    Qtilde = generateQtilde(pathinfo,MPC_vars,ModelParams,Xk,i);
    % add omega regularization
    if i == MPC_vars.N+1
        Qtilde(ModelParams.stateindex_omega,ModelParams.stateindex_omega) = MPC_vars.qOmegaNmult*MPC_vars.qOmega;
    else
        Qtilde(ModelParams.stateindex_omega,ModelParams.stateindex_omega) = MPC_vars.qOmega;
    end
    % make Qtilde symetric (not symetric due to numerical issues)
    Qtilde = 0.5 *(Qtilde+Qtilde');
    % Qk = contouring-lag error and real-input cost
    Qk = 2*blkdiag(Qtilde,diag([MPC_vars.rD,MPC_vars.rDelta,MPC_vars.rVtheta]));  %%%%%diag([MPC_vars.rD,MPC_vars.rDelta,MPC_vars.rVtheta])为控制量的权重
    % scale cost
    Qk = blkdiag(MPC_vars.invTx,MPC_vars.invTu)*Qk*blkdiag(MPC_vars.invTx,MPC_vars.invTu) + 1e-12*eye(10);  %%矩阵归一化，提高计算效率
end

% compute linear contouring and lag errors
function Qtilde = generateQtilde(pathinfo,MPC_vars,ModelParams,Xk,i)
    if i == MPC_vars.N+1
        Q = diag([MPC_vars.qCNmult*MPC_vars.qC, MPC_vars.qL]);            %[terminal contouring cost    lag cost] qCNmult为contoring cost权重调节系数，预测的最后一步加大权重，如何理解呢？
    else
        Q = diag([MPC_vars.qC, MPC_vars.qL]);
    end
        
    theta_virt=mod(Xk(end),pathinfo.ppx.breaks(end));                     %与最大值取余是为了循环下去
    [grad_eC, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams,Xk(1), Xk(2));
    errorgrad = [grad_eC; grad_eL]; 
    Qtilde = errorgrad'*Q*errorgrad; 
end

function [grad_eC, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams, x_phys, y_phys)

    [deC_dtheta, deL_dtheta, cos_phi_virt, sin_phi_virt] = getderror_dtheta(pathinfo, theta_virt, x_phys, y_phys);
    
    grad_eC = [ sin_phi_virt, -cos_phi_virt, zeros(1, ModelParams.nx-3), deC_dtheta];
    grad_eL = [-cos_phi_virt, -sin_phi_virt, zeros(1, ModelParams.nx-3), deL_dtheta];
end

function [deC_dtheta, deL_dtheta, cos_phi_virt, sin_phi_virt] = getderror_dtheta(pathinfo, theta_virt, x_phys, y_phys)
    dxvirt_dtheta=ppval(pathinfo.dppx,theta_virt); %d x_virt / d theta
    dyvirt_dtheta=ppval(pathinfo.dppy,theta_virt); %d y_virt / d theta
    
    phi_virt=atan2(dyvirt_dtheta,dxvirt_dtheta); %orientation of virtual position
    % virtual positions
    x_virt=ppval(pathinfo.ppx,theta_virt);
    y_virt=ppval(pathinfo.ppy,theta_virt);
    
    % difference in position between virtual and physical
    Dx=x_phys-x_virt;
    Dy=y_phys-y_virt;

    dphivirt_dtheta=getdphivirt_dtheta(theta_virt,pathinfo);

    cos_phi_virt=cos(phi_virt);
    sin_phi_virt=sin(phi_virt);

    tmp1=[dphivirt_dtheta, 1];
    tmp2=[cos_phi_virt; sin_phi_virt];

    MC=[ Dx Dy; dyvirt_dtheta -dxvirt_dtheta];
    ML=[-Dy Dx; dxvirt_dtheta  dyvirt_dtheta];

    deC_dtheta = tmp1 * MC * tmp2;
    deL_dtheta = tmp1 * ML * tmp2;
end

%%%%获得偏航角对曲线sita的导数
function dphivirt_dtheta=getdphivirt_dtheta(theta_virt,pathinfo)
    % computes {d phi_virt / d theta} evaluated at theta_k

    dxdth=ppval(pathinfo.dppx,theta_virt); %d x_virt / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); %d y_virt / d theta
    d2xdth2=ppval(pathinfo.ddppx,theta_virt); %d2 x_virt / d theta2
    d2ydth2=ppval(pathinfo.ddppy,theta_virt); %d2 y_virt / d theta2

    numer=dxdth*d2ydth2 - dydth*d2xdth2;
    denom=dxdth^2 + dydth^2;

    dphivirt_dtheta=numer/denom;
end

% GENERATING f
function f = generatef(pathinfo,MPC_vars,ModelParams,Xk,i)

    x_phys = Xk(1);
    y_phys = Xk(2);

    theta_virt=mod(Xk(end),pathinfo.ppx.breaks(end));
    [eC, eL] = getErrors(pathinfo, theta_virt,x_phys,y_phys);
    e=[eC;eL];
    [grad_eC, grad_eL] = getErrorGradient(pathinfo, theta_virt, ModelParams, x_phys, y_phys);
    grad_e = [grad_eC; grad_eL];
    
    if i == MPC_vars.N+1
        Q = diag([MPC_vars.qCNmult*MPC_vars.qC, MPC_vars.qL]);
    else
        Q = diag([MPC_vars.qC, MPC_vars.qL]);                             %轨迹误差法线和切线的权重
    end
  
    fx=2*e'*Q*grad_e - 2*Xk'*grad_e'*Q*grad_e;
    fT = [fx, zeros(1,ModelParams.nu-1), -MPC_vars.qVtheta];
    f=fT';
    
    f = blkdiag(MPC_vars.invTx,MPC_vars.invTu)*f;
end

function [eC, eL] = getErrors(pathinfo, theta_virt,x_phys,y_phys)
    dxdth=ppval(pathinfo.dppx,theta_virt); % d x / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); % d y / d theta

    % virtual positions
    x_virt=ppval(pathinfo.ppx,theta_virt); 
    y_virt=ppval(pathinfo.ppy,theta_virt);
    
    phi_virt=atan2(dydth,dxdth);
    
    % define these to reduce calls to trig functions
    sin_phi_virt = sin(phi_virt);
    cos_phi_virt = cos(phi_virt);

    % contouring and lag error estimates
    eC = -sin_phi_virt*(x_virt - x_phys) + cos_phi_virt*(y_virt - y_phys); %法线和切线的偏差from 公式（11a and 11b）
    eL =  cos_phi_virt*(x_virt - x_phys) + sin_phi_virt*(y_virt - y_phys);
   
end

% EQUALITY CONSTRAINTS  动力学方程
function [Ak,Bk,gk] = getEqualityConstraints(Xk,Uk,MPC_vars,ModelParams) 

    nx = ModelParams.nx;
    nu = ModelParams.nu;
    % linearize and discretize nonlinear bicycle model
    [Ad, Bd, gd]=DiscretizedLinearizedModel(Xk,Uk,ModelParams,MPC_vars.Ts);
    % constructing augmented system with state-input scaling
    % Xk+1=Ak*Xk+Bk*det_Uk+gk
    Ak = [MPC_vars.Tx*Ad*MPC_vars.invTx MPC_vars.Tx*Bd*MPC_vars.invTu; zeros(nu,nx) eye(nu)];
    Bk = [MPC_vars.Tx*Bd*MPC_vars.invTu;eye(nu)];
    gk = [MPC_vars.Tx*gd;zeros(nu,1)];
    
end

% INEQUALITY CONSTRAINTS
function [Ck, ug, lg] = getInequalityConstraints(border,MPC_vars,ModelParams)  %三角形不等式
    nx = ModelParams.nx;
    nu = ModelParams.nu;
    x1=border(1);
    y1=border(2);
    x2=border(3);
    y2=border(4);
    
    % numerator and denominator of slope of border. m = - (x2-x1)/(y2-y1)
    numer=-(x2-x1);
    denom=(y2-y1);

    dbmax=max(numer*x1-denom*y1,numer*x2-denom*y2);
    dbmin=min(numer*x1-denom*y1,numer*x2-denom*y2);

    Ck=zeros(1,nx+nu);
    Ck(1,1:2)=[numer -denom];
%     Ck(1,7)=1;
    ug= dbmax;
    lg= dbmin;
    
    Ck = Ck*blkdiag(MPC_vars.invTx,MPC_vars.invTu);
    
    
end
%%
%ObstaclesAvodance
function [Dk, dk] = getInequalityConstraints_OA(MPC_vars,ModelParams,grid_isOccupiedBy,X_all,SimN)  %椭圆不等式
    nx = ModelParams.nx;
    nu = ModelParams.nu;
    a=ModelParams.L *1.3;
    b=ModelParams.W*1.3; 
   
    Dk=zeros(4,nx+nu); 
    dk=zeros(4,1);
    if max(max(abs(grid_isOccupiedBy)))==0 %in C -1, first car is car0    %如果周边没车，那么就默认路面边界 
        % if all elements of grid_isOccupiedBy are zeros (ie. no cars nearby),
        % then use default borders
        
    else                                                                  %否则，重新规划边界 
%   
     
        for i= 1:4
            det_X=X_all(1,1)-X_all(1,i+1);
            det_Y=X_all(2,1)-X_all(2,i+1);
            det_T=[det_X;det_Y];
            phi=X_all(3,i+1);
            A_T=[b*cos(phi) b*sin(phi);-a*sin(phi) a*cos(phi)]/a/b;
            d_A_T=[-b*sin(phi) b*cos(phi);-a*cos(phi) -a*sin(phi)]/a/b;
            drou_dxANDdy=2*A_T'*A_T*det_T;
            drou_dphi=2*(A_T*det_T)'*(d_A_T*det_T);
            Dk(i,1:3)=-[drou_dxANDdy' drou_dphi*0];
            rou0=((A_T(1,:)*det_T)^2+(A_T(2,:)*det_T)^2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dk(i,1)=rou0-1+Dk(i,1:3)*X_all(1:3,1)*0;      
            
        end

    end

    if SimN <5
       dk=dk*0; 
       Dk=Dk*blkdiag(MPC_vars.invTx,MPC_vars.invTu)*0;  
    end
end

function [Ek, ek] =  getInequalityConstraints_OB(pathinfo,MPC_vars,ModelParams,Xk) %v^2/R<0.4g约束 过弯约束 
   nx = ModelParams.nx;
    nu = ModelParams.nu;

    theta_virt=mod(Xk(end),pathinfo.ppx.breaks(end));
    dxdth=ppval(pathinfo.dppx,theta_virt); % d x / d theta
    dydth=ppval(pathinfo.dppy,theta_virt); % d y / d theta
    ddxdth=ppval(pathinfo.ddppx,theta_virt); % d x / d theta
    ddydth=ppval(pathinfo.ddppy,theta_virt); % d y / d theta
    rou =(dxdth*ddydth-ddxdth*dydth)/(dxdth^2+dydth^2)^(3/2);
    Ek=zeros(1,nx+nu); 
    ek=0;
    
    if rou > 1/100
        R=1/rou;
        Ek(:,4)=1;
        ek=sqrt(0.35*9.81*R);
        Ek = Ek*blkdiag(MPC_vars.invTx,MPC_vars.invTu);
    end
end

% BOUNDS
function [lb, ub]=getBounds(MPC_vars,ModelParams)

lb = MPC_vars.bounds(:,1);
ub = MPC_vars.bounds(:,2);


end
