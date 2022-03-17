
%% MPCC Simulation Script ģ��Ԥ���������� ����Ŀ���ǵ�ǰ�켣��Ŀ��켣����/���߾���ֵ֮����С���滮�켣��Ŀ��켣�ϵ�ͶӰ��󡢿�����������С
clear
close all 
clc 

kk=1/310;  %kk=0Ϊֱ����ʻ��kk=1Ϊcircle
%% add spline library
addpath('splines');   
% addpath('~/Documents/GitHub/hpipm/interfaces/matlab/hpipm_matlab')
%% Load Parameters
% CarModel = 'ORCA';
CarModel = 'FullSize';

MPC_vars = getMPC_vars(CarModel);                                         %MPC����
ModelParams=getModelParams(MPC_vars.ModelNo);                             %��������
% choose optimization interface options: 'Yalmip','CVX','hpipm','quadprog'
MPC_vars.interface = 'quadprog';
 
nx = ModelParams.nx;                              %%%%%%%%%%%%%%%%%%%%%%%%%״̬��������[ X Y yawAngle vx vy omega virtual_position]   virtual position Ԥ��켣�ĵ�����Ŀ��켣��ͶӰ
nu = ModelParams.nu;                              %%%%%%%%%%%%%%%%%%%%%%%%%���Ʊ�������[ ��������ź�  ������ת�� ����]
N = MPC_vars.N;                                   %%%%%%%%%%%%%%%%%%%%%%%%%Ԥ�ⲽ��
Ts = MPC_vars.Ts;                                 %%%%%%%%%%%%%%%%%%%%%%%%%���㲽�� 
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
plot(track.outer(1,:),track.outer(2,:),'r','linewidth',3)                                 %��ɫΪԤ����������*0.75���İ�ȫ��·�ߣ�Ҳ������������Ĺ�ʽ��6d����Լ��
hold on 
plot(track.inner(1,:),track.inner(2,:),'r','linewidth',3)                               %��ɫΪԤ����������*0.75���İ�ȫ��·��
% plot(track2.outer(1,:),track2.outer(2,:),'k')                             %ʵ�ʵ�·��
% plot(track2.inner(1,:),track2.inner(2,:),'k','linewidth',3)                             %ʵ�ʵ�·�� 
%    axis([0,100,25,50,]);
%% Simulation lenght and plotting 
simN = 150;                                                         %�����ܲ���           
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
[traj, borders] =splinify(track);                                         %���ο��켣�����ߺͱ߽��ߣ����ߣ�������ϵ������һ�׶���΢��ϵ��
tl = traj.ppy.breaks(end);                                                 %·���ܳ�    
  
% store all data in one struct
TrackMPC = struct('traj',traj,'borders',borders,'track_center',track.center,'tl',tl);
%% Define starting position ���ҵ���ʼλ�õ�
startIdx = 310*kk+1+1; %point (in terms of track centerline array) allong the track  %$##$$#%%%^#!@#$%^&*()@#$%^&*(!@#$%^&*()@#$%^&*(@#$%^&*()@#$%^&*(@#$%^&*)(*&^%$#@(*&^%$#@@#$%^&*()PLMNGFEW#$%^&*&^%$#@!@#$%^&*(*&^%$#@!@#$%^&*(*&^%$#@
% where the car starts, on the center line, aligned with the track, driving
% straight with vx0
%since the used bicycle model is not well defined for slow velocities use vx0 > 0.5
if strcmp(CarModel,'ORCA')                                                %ѡ��ģ�Ͷ�Ӧ�ĳ�ʼ����
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
%% First initial guess Ȼ��ѳ�ʼ��NP��������������������ʼ��
x = repmat(x0,1,N+1); % all points identical to current measurment

% first inital guess, all points on centerline aligned with centerline
% spaced as the car would drive with vx0
for i = 2:N+1
    theta_next = x(ModelParams.stateindex_theta,i-1)+Ts*vx0;              %�������������·�̵㣻 ����'ModelParams.stateindex_theta'Ϊ����x�ĵڼ��У�Ҳ���ǵڼ���״̬����x�ĳ�ʼֵ���������x0
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
plot(track.outer(1,:),track.outer(2,:),'r')                               %��ɫΪԤ����һ������İ�ȫ��·�� 
hold on
plot(track.inner(1,:),track.inner(2,:),'r')                               %��ɫΪԤ����һ������İ�ȫ��·��
plot(track2.outer(1,:),track2.outer(2,:),'k')                             %ʵ�ʵ�·��
plot(track2.inner(1,:),track2.inner(2,:),'k')                             %ʵ�ʵ�·��  
u = zeros(3,N); % zero inputs
uprev = zeros(3,1); % last input is zero 
%% Ohter cars ����ϰ���
Y0 = zeros(nx,n_cars-1);
closestID_last0=[];
kkk0=zeros(1,4);
[Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,0,Y0,closestID_last0,kkk0);                          %ÿһ���ϰ�����X״̬ʸ��
kkk0=kkk;
if size(Y,2) ~= n_cars-1
    error('n_cars and the number of obstacles in "Y" does not match')
end 
%% Initialize logging arrays                                     
X_log = zeros(nx*(N+1),simN); 
U_log = zeros(nu*N,simN);
B_log = zeros(4*N,simN);                                                  %�߽���־
qpTime_log = zeros(1,simN); 
%% initializtion
% solve problem 5 times without applying input
% inspiered by sequential quadratic programming (SQP) 
convergeflag=1;
for i = 1:5   
%    formulate MPCC problem and solve it 
    Iter_damping = 0.5; % 0 no damping
    [x_up, u_up, b, exitflag,info] = optimizer_self(TrackMPC,MPC_vars,ModelParams, n_cars, Y, x, u, x0, uprev,i, convergeflag);
    [Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,i,Y0,closestID_last0,kkk0);                          %ÿһ���ϰ�����X״̬ʸ��
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
%% Simulation ������ 
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
    [Y, closestID_last,kkk]= ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,i,Y0,closestID_last0,kkk0);                          %ÿһ���ϰ�����X״̬ʸ��
    kkk0=kkk;    
    x0 = SimTimeStep(x(:,1),u(:,1),Ts,ModelParams)';                      %���±���ϵͳ״̬��
    x0 = unWrapX0(x0);                                                    %ƫ���Ƿ�Χ������-pi pi��
    [ theta, last_closestIdx] = findTheta(x0,track.center,traj.ppx.breaks,trackWidth,last_closestIdx);
    x0(ModelParams.stateindex_theta) = theta ;
    uprev = u(:,1);   
%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% plotting and logging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if plotOn == 1  
        PlotPrediction(x,u,b,Y,track,track2,traj,MPC_vars,ModelParams,i,qpTime_log(i))     %xΪԤ��Ĺ켣��uΪ��������bΪ·��boundary��YΪ�����ϰ�����
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
    LapTime = diff(LapStart)*Ts;                                           %��һȦ����Ҫ��ʱ��                                 
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