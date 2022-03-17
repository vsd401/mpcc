
%% MPCC Simulation Script 模型预测轮廓控制 控制目标是当前轨迹和目标轨迹切线/法线绝对值之差最小、规划轨迹在目标轨迹上的投影最大、控制量消耗最小
clear
close all 
clc 

kk=1/310;  %kk=0为直线行驶；kk=1为circle
%% add spline library
addpath('splines');   
% addpath('~/Documents/GitHub/hpipm/interfaces/matlab/hpipm_matlab')
%% Load Parameters
% CarModel = 'ORCA';
CarModel = 'FullSize';

MPC_vars = getMPC_vars(CarModel);                                         %MPC参数
ModelParams=getModelParams(MPC_vars.ModelNo);                             %整车参数
% choose optimization interface options: 'Yalmip','CVX','hpipm','quadprog'
MPC_vars.interface = 'quadprog';
 
nx = ModelParams.nx;                              %%%%%%%%%%%%%%%%%%%%%%%%%状态变量个数[ X Y yawAngle vx vy omega virtual_position]   virtual position 预测轨迹的导数在目标轨迹的投影
nu = ModelParams.nu;                              %%%%%%%%%%%%%%%%%%%%%%%%%控制变量个数[ 电机控制信号  方向盘转角 车速]
N = MPC_vars.N;                                   %%%%%%%%%%%%%%%%%%%%%%%%%预测步长
Ts = MPC_vars.Ts;                                 %%%%%%%%%%%%%%%%%%%%%%%%%计算步长 
%% import an plot track 
% use normal ORCA Track
% % % % % % % % % % % load Tracks/track2.mat
% use RCP track  
% % % % % % % % % % % % % % % % load Tracks/trackMobil.mat 
% % % % % % % % % % % % % % % % track2 = trackMobil;
load Tracks/Trackdata.mat 
track2 = Trackdata;
% shrink track by half of the car widht plus safety margin  
% TODO implement orientation depending shrinking in the MPC constraints
safteyScaling = 1.5;  
[track,track2] = borderAdjustment(track2,ModelParams,safteyScaling);
track.outer(1,:)=track.outer(1,:)+50*kk; 
track.inner(1,:)=track.inner(1,:)+50*kk; 
track.center(1,:)=track.center(1,:)+50*kk;
track2.outer(1,:)=track2.outer(1,:)+50*kk;
track2.inner(1,:)=track2.inner(1,:)+50*kk;
track2.center(1,:)=track2.center(1,:)+50*kk; 

trackWidth = norm(track.inner(:,1)-track.outer(:,1));   
% plot shrinked and not shrinked track   
figure(1);   
plot(track.outer(1,:),track.outer(2,:),'r','linewidth',3)                                 %红色为预留出（车宽*0.75）的安全道路线，也就是论文里面的公式（6d）的约束
hold on 
plot(track.inner(1,:),track.inner(2,:),'r','linewidth',3)                               %红色为预留出（车宽*0.75）的安全道路线
% plot(track2.outer(1,:),track2.outer(2,:),'k')                             %实际道路线
% plot(track2.inner(1,:),track2.inner(2,:),'k','linewidth',3)                             %实际道路线 
%    axis([0,100,25,50,]);
%% Simulation lenght and plotting 
simN = 150;                                                         %计算总步数           
%0=no plots, 1=plot predictions
plotOn = 1; 
%0=real time iteration, 1=fixed number of QP iterations, 2=fixed number of damped QP iterations
QP_iter = 2;
% number of cars  
% there are two examples one with no other cars and one with 4 other cars
% (inspired by the set up shown in the paper)
% n_cars = 1; % no other car 
n_cars = 5; % 4 other cars
%% Fit spline to track  
% TODO spline function only works with regular spaced points.
% Fix add function which given any center line and bound generates equlally
% space tracks.
[traj, borders] =splinify(track);                                         %求解参考轨迹中心线和边界线（红线）的样条系数及其一阶二阶微分系数
tl = traj.ppy.breaks(end);                                                 %路程总长    
  
% store all data in one struct
TrackMPC = struct('traj',traj,'borders',borders,'track_center',track.center,'tl',tl);
%% Define starting position 先找到起始位置点
startIdx = 310*kk+1+1; %point (in terms of track centerline array) allong the track  %$##$$#%%%^#!@#$%^&*()@#$%^&*(!@#$%^&*()@#$%^&*(@#$%^&*()@#$%^&*(@#$%^&*)(*&^%$#@(*&^%$#@@#$%^&*()PLMNGFEW#$%^&*&^%$#@!@#$%^&*(*&^%$#@!@#$%^&*(*&^%$#@
% where the car starts, on the center line, aligned with the track, driving
% straight with vx0
%since the used bicycle model is not well defined for slow velocities use vx0 > 0.5
if strcmp(CarModel,'ORCA')                                                %选择模型对应的初始车速
    vx0 = 1;
elseif strcmp(CarModel,'FullSize') 
    vx0 = 15;
end   
   
% find theta that coresponds to the 10th point on the centerline
[theta, ~] = findTheta([track.center(1,startIdx),track.center(2,startIdx)],track.center,traj.ppx.breaks,trackWidth,startIdx);

x0 = [track.center(1,startIdx),track.center(2,startIdx),... % point on centerline
      atan2(ppval(traj.dppy,theta),ppval(traj.dppx,theta)),... % aligned with centerline
      vx0 ,0,0,theta]'; %driving straight with vx0, and correct theta progress
    
% the find theta function performs a local search to find the projection of
% the position onto the centerline, therefore, we use the start index as an
% starting point for this local search
last_closestIdx = startIdx;
%% First initial guess 然后把初始的NP步迭代出来，还是做初始化
x = repmat(x0,1,N+1); % all points identical to current measurment

% first inital guess, all points on centerline aligned with centerline
% spaced as the car would drive with vx0
for i = 2:N+1
    theta_next = x(ModelParams.stateindex_theta,i-1)+Ts*vx0;              %该命令行是求解路程点； 其中'ModelParams.stateindex_theta'为表述x的第几行，也就是第几个状态量，x的初始值来自上面的x0
    phi_next = atan2(ppval(traj.dppy,theta_next),ppval(traj.dppx,theta_next));%
    % phi_next can jump by two pi, make sure there are no jumps in the
    % initial guess
    if (x(ModelParams.stateindex_phi,i-1)-phi_next) < -pi
        phi_next = phi_next-2*pi;
    end
    if (x(ModelParams.stateindex_phi,i-1)-phi_next) > pi
        phi_next = phi_next+2*pi;
    end 
    x(:,i) = [ppval(traj.ppx,theta_next),ppval(traj.ppy,theta_next),... % point on centerline
              phi_next,... % aligned with centerline
              vx0 ,0,0,theta_next]'; %driving straight with vx0, and correct theta progress
end
figure(77)
plot(x(1,:),x(2,:),'*')   
hold on
plot(track.outer(1,:),track.outer(2,:),'r')                               %红色为预留出一个车宽的安全道路线 
hold on
plot(track.inner(1,:),track.inner(2,:),'r')                               %红色为预留出一个车宽的安全道路线
plot(track2.outer(1,:),track2.outer(2,:),'k')                             %实际道路线
plot(track2.inner(1,:),track2.inner(2,:),'k')                             %实际道路线  
u = zeros(3,N); % zero inputs
uprev = zeros(3,1); % last input is zero 
%% Ohter cars 添加障碍车
Y0 = zeros(nx,n_cars-1);
closestID_last0=[];
kkk0=zeros(1,4);
[Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,0,Y0,closestID_last0,kkk0);                          %每一个障碍车的X状态矢量
kkk0=kkk;
if size(Y,2) ~= n_cars-1
    error('n_cars and the number of obstacles in "Y" does not match')
end 
%% Initialize logging arrays                                     
X_log = zeros(nx*(N+1),simN); 
U_log = zeros(nu*N,simN);
B_log = zeros(4*N,simN);                                                  %边界日志
qpTime_log = zeros(1,simN); 
%% initializtion
% solve problem 5 times without applying input
% inspiered by sequential quadratic programming (SQP) 
convergeflag=1;
for i = 1:5   
%    formulate MPCC problem and solve it 
    Iter_damping = 0.5; % 0 no damping
    [x_up, u_up, b, exitflag,info] = optimizer_self(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev,i, convergeflag);
    [Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,i,Y0,closestID_last0,kkk0);                          %每一个障碍车的X状态矢量
    kkk0=kkk;
    Y0=Y;
    closestID_last0=closestID_last;
    x = Iter_damping*x + (1-Iter_damping)*x_up;
    u = Iter_damping*u + (1-Iter_damping)*u_up;

        if plotOn == 1      
    %         plot predictions 
            if   i > 0
                 PlotPrediction(x,u,b,Y,track,track2,traj,MPC_vars,ModelParams,i,info.QPtime)   
            end
        end    
end   
%% Simulation 主程序 
for i = 2: simN   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%% MPCC-Call %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % augment state and inputs by shifting previus optimal solution 
    [x,u] = augState(x,u,x0,MPC_vars,ModelParams,tl);   
    %  formulate MPCC problem and solve it  
    if QP_iter == 0   
        [x, u, b, exitflag,info] = optimizer(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev);
        
        qpTime_log(i) = info.QPtime; 
    elseif QP_iter == 1     
        % doing multiple "SQP" steps
        for k = 1:2
            [x, u, b, exitflag,info] = optimizer(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev);
            qpTime_log(i) = qpTime_log(i) + info.QPtime;
        end
    elseif QP_iter == 2  
        % doing multiple damped "SQP" steps
        for k = 1:2 
            Iter_damping = 0.75; % 0 no damping  
            [x_up, u_up, b, exitflag,info] = optimizer_self(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev,i, convergeflag);
            convergeflag=exitflag;
            x = Iter_damping*x + (1-Iter_damping)*x_up;    
            u = Iter_damping*u + (1-Iter_damping)*u_up;  
            qpTime_log(i) = qpTime_log(i) + info.QPtime;   
        end 
    else 
        error('invalid QP_iter value')
    end
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% simulate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    [Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,i,Y0,closestID_last0,kkk0);                          %每一个障碍车的X状态矢量
    kkk0=kkk;    
    x0 = SimTimeStep(x(:,1),u(:,1),Ts,ModelParams)';                      %更新被控系统状态量
    x0 = unWrapX0(x0);                                                    %偏航角范围纠正【-pi pi】
    [ theta, last_closestIdx] = findTheta(x0,track.center,traj.ppx.breaks,trackWidth,last_closestIdx);
    x0(ModelParams.stateindex_theta) = theta ;
    uprev = u(:,1);   
%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% plotting and logging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if plotOn == 1  
        PlotPrediction(x,u,b,Y,track,track2,traj,MPC_vars,ModelParams,i,qpTime_log(i))     %x为预测的轨迹；u为控制量；b为路的boundary；Y为其他障碍车；
    end
    
    % log predictions and time 
    X_log(:,i) = reshape(x,(N+1)*7,1);
    U_log(:,i) = reshape(u,(N)*3,1);
    B_log(:,i) = reshape(b,N*4,1); 
    if n_cars>1
        for j=1:n_cars-1
            Y(:,j)=unWrapX0(Y(:,j));
        end
        Y0=unWrapX0(Y);
    end 
    closestID_last0=closestID_last;         
    
    
end  


PlotLog( X_log,U_log,Y,track,track2,simN,Ts) 

%% Generating Stats
a = 1;
LapStart= zeros(1,simN);
for i=1:simN-1
    if X_log(ModelParams.stateindex_theta,i+1) - X_log(ModelParams.stateindex_theta,i) < -0.9*tl
        LapStart(a) = i;
        a = a+1;
    end
end

if length(LapStart) > 1
    LapTime = diff(LapStart)*Ts;                                           %跑一圈所需要的时间                                 
else
    LapTime = NaN;
end

disp('------------------------------------')
disp(['Lap Time(s): ',num2str(LapTime)])
disp('------------------------------------')
disp(['Mean Computation Time: ',num2str(mean(qpTime_log))])
disp(['Max Computation Time: ',num2str(max(qpTime_log))])
disp(['Min Computation Time: ',num2str(min(qpTime_log))])
disp('------------------------------------')