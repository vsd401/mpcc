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

function [Y, closestID_last,kkk] = ObstacelsState(track,traj,trackWidth,n_cars,kk,MPC_vars,SimN,Y0,closestID_last0,kkk0)

if n_cars == 1 
    Y = [];
    closestID_last=[];
    kkk=zeros(1,4);
elseif n_cars == 5
    
    
    
    initalV = [15, 15, 15, 15]*0.1;   %障碍车初始车速
    for i=n_cars-1
        vx0_Op1=initalV(1);
        vx0_Op2=initalV(2);
        vx0_Op3=initalV(3);
        vx0_Op4=initalV(4);
    end
    
    
     if SimN <3
    % find theta that coresponds to the 10th point on the centerline
        startIdx_Op1 = 10+310*kk;  
        [theta_Op1, closestIdx1] = findTheta([track.center(1,startIdx_Op1),track.center(2,startIdx_Op1)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op1); %
        kkk(1)=0;
    else
        [theta_Op1, closestIdx1] = findTheta([track.center(1,closestID_last0(1))+vx0_Op1*MPC_vars.Ts*cos(Y0(3,1)),track.center(2,closestID_last0(1))+vx0_Op1*MPC_vars.Ts*sin(Y0(3,1))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(1));
        if closestIdx1  >= closestID_last0(1)  && vx0_Op1>=0
            kkk(1)=kkk0(1)+1;
            [theta_Op1, closestIdx1] = findTheta([track.center(1,closestID_last0(1))+kkk(1)*vx0_Op1*MPC_vars.Ts*cos(Y0(3,1)),track.center(2,closestID_last0(1))+kkk(1)*vx0_Op1*MPC_vars.Ts*sin(Y0(3,1))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(1));
            if closestIdx1 > closestID_last0(1)
               kkk(1)=0;
            end
        end 
     end 
    
    x0_Op1 = [track.center(1,closestIdx1)+(-trackWidth/3)*sin(-Y0(3,1)),track.center(2,closestIdx1)+(-trackWidth/3)*cos(-Y0(3,1)),... % point on centerline       [X  Y+半个轮距  phi vx vy d_phi theta_Op1]
          atan2(ppval(traj.dppy,theta_Op1),ppval(traj.dppx,theta_Op1)),... % aligned with centerline   
          vx0_Op1 ,0,0,theta_Op1]'; %driving straight with vx0, and correct theta progress

      
    if SimN <3
         %     vx0_Op2 = 2; 
         % find theta that coresponds to the 10th point on the centerline
        startIdx_Op2 = 20+310*kk; 
        [theta_Op2,closestIdx2] = findTheta([track.center(1,startIdx_Op2),track.center(2,startIdx_Op2)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op2);
        kkk(2)=0;
    else
        [theta_Op2, closestIdx2] = findTheta([track.center(1,closestID_last0(2))+vx0_Op1*MPC_vars.Ts*cos(Y0(3,2)),track.center(2,closestID_last0(2))+vx0_Op1*MPC_vars.Ts*sin(Y0(3,2))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(2));
        if closestIdx2  >= closestID_last0(2) && vx0_Op1>=0
            kkk(2)=kkk0(2)+1;
            [theta_Op2, closestIdx2] = findTheta([track.center(1,closestID_last0(2))+kkk(2)*vx0_Op1*MPC_vars.Ts*cos(Y0(3,2)),track.center(2,closestID_last0(2))+kkk(2)*vx0_Op1*MPC_vars.Ts*sin(Y0(3,2))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(2));
            if closestIdx2 > closestID_last0(2)
               kkk(2)=0;
            end
        end
    end  
    x0_Op2 = [track.center(1,closestIdx2)+(trackWidth/3)*sin(-Y0(3,2)),track.center(2,closestIdx2)+(trackWidth/3)*cos(-Y0(3,2)),... % point on centerline
          atan2(ppval(traj.dppy,theta_Op2),ppval(traj.dppx,theta_Op2)),... % aligned with centerline
          vx0_Op2 ,0,0,theta_Op2]'; %driving straight with vx0, and correct theta progress

      
    if SimN <3
         startIdx_Op3 = 33+310*kk; 
         kkk(3)=0;
        %     vx0_Op3 = 2; 
        % find theta that coresponds to the 10th point on the centerline
        [theta_Op3, closestIdx3] = findTheta([track.center(1,startIdx_Op3),track.center(2,startIdx_Op3)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op3);
    else
        [theta_Op3, closestIdx3] = findTheta([track.center(1,closestID_last0(3))+vx0_Op1*MPC_vars.Ts*cos(Y0(3,3)),track.center(2,closestID_last0(3))+vx0_Op1*MPC_vars.Ts*sin(Y0(3,3))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(3));
         if closestIdx3  >= closestID_last0(3) && vx0_Op1>=0
            kkk(3)=kkk0(3)+1;
            [theta_Op3, closestIdx3] = findTheta([track.center(1,closestID_last0(3))+kkk(3)*vx0_Op1*MPC_vars.Ts*cos(Y0(3,3)),track.center(2,closestID_last0(3))+kkk(3)*vx0_Op1*MPC_vars.Ts*sin(Y0(3,3))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(3));
             if closestIdx3 > closestID_last0(3)
                kkk(3)=0;
            end
         end
            
    end     
    x0_Op3 = [track.center(1,closestIdx3)+(trackWidth/3)*sin(-Y0(3,3)),track.center(2,closestIdx3)+(trackWidth/3)*cos(-Y0(3,3)),... % point on centerline
          atan2(ppval(traj.dppy,theta_Op3),ppval(traj.dppx,theta_Op3)),... % aligned with centerline
          vx0_Op3 ,0,0,theta_Op3]'; %driving straight with vx0, and correct theta progress

    
    if SimN <3
        startIdx_Op4 = 31+310*kk; 
    %     vx0_Op4 = 2; 
    % find theta that coresponds to the 10th point on the centerline
        [theta_Op4, closestIdx4] = findTheta([track.center(1,startIdx_Op4),track.center(2,startIdx_Op4)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op4);
        kkk(4)=0;
    else
        [theta_Op4, closestIdx4] = findTheta([track.center(1,closestID_last0(4))+vx0_Op1*MPC_vars.Ts*cos(Y0(3,4)),track.center(2,closestID_last0(4))+vx0_Op1*MPC_vars.Ts*sin(Y0(3,4))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(4));
        if closestIdx4  >= closestID_last0(4) && vx0_Op1>0
            kkk(4)=kkk0(4)+1;
            [theta_Op4, closestIdx4] = findTheta([track.center(1,closestID_last0(4))+kkk(4)*vx0_Op1*MPC_vars.Ts*cos(Y0(3,4)),track.center(2,closestID_last0(4))+kkk(4)*vx0_Op1*MPC_vars.Ts*sin(Y0(3,4))],track.center,traj.ppx.breaks,trackWidth,closestID_last0(4));
            if closestIdx4 > closestID_last0(4)
               kkk(4)=0;
            end
        end
    end    
    x0_Op4 = [track.center(1,closestIdx4)+(-trackWidth/4*0)*sin(-Y0(3,4)),track.center(2,closestIdx4)+(-trackWidth/4*0)*cos(-Y0(3,4)),... % point on centerline
          atan2(ppval(traj.dppy,theta_Op4),ppval(traj.dppx,theta_Op4)),... % aligned with centerline
          vx0_Op4 ,0,0,theta_Op4]'; %driving straight with vx0, and correct theta progress
    x0_Op1 = unWrapX0(x0_Op1); 
    x0_Op2 = unWrapX0(x0_Op2); 
    x0_Op3 = unWrapX0(x0_Op3); 
    x0_Op4 = unWrapX0(x0_Op4); 
    Y = [x0_Op1,x0_Op2,x0_Op3,x0_Op4];
    closestID_last=[closestIdx1 closestIdx2 closestIdx3 closestIdx4 ]; 
else
    error('only 0 or 4 obstacles are pre programmed settings, feel free to change the obstacle constellation')
end