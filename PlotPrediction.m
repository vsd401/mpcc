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

function [  ] = PlotPrediction( x,u,b,Y,track,track2,traj,MPC_vars,ModelParams,simN,simTime )
    
    N = MPC_vars.N; 
    Ts = MPC_vars.Ts;
    
    tl = traj.ppx.breaks(end);
    
    figure(1);
    width=800;%宽度，像素数
    height=270;%高度
    left=0;%距屏幕左下角水平距离
    bottem=200;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height],'color','w')
    plot(track.outer(1,:),track.outer(2,:),'r','linewidth',3)
%     axis([-50,50,25,50,]);
    hold on
    plot(track.inner(1,:),track.inner(2,:),'r','linewidth',3)
%     plot(track2.outer(1,:),track2.outer(2,:),'k')
%     plot(track2.inner(1,:),track2.inner(2,:),'k')
%     plot(b(:,1),b(:,2),'black.-')
%     plot(b(:,3),b(:,4),'black.-')
    plot(ppval(traj.ppx,mod(x(7,:),tl)),ppval(traj.ppy,mod(x(7,:),tl)),':k')
    plot(x(1,:),x(2,:),'--')
    carBox(x(:,1),ModelParams.W/2,ModelParams.L/2)
    if ~isempty(Y)
        for i=1:size(Y,2)
            carBox(Y(:,i),ModelParams.W/2,ModelParams.L/2)
            aa=ModelParams.L *1.1;
            bb=ModelParams.W*1.1; 
               z=linspace(0,2*pi,100);
               xx=aa*sin(z);
               yy=bb*cos(z);
               phi=Y(3,i);
               A=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[xx;yy]+kron(ones(size(z)),Y(1:2,i));
               plot(A(1,:),A(2,:))                            %画椭圆
        end
    end
    xlabel('X [m]')
    ylabel('Y [m]')
   axis([0+x(7,1)*0,100*0+x(7,1)*0+73,-35,10]); 
%       axis([0+x(7,1)-3,80+x(7,1),28,48,]);   
%             plot(X_all(1,1),X_all(2,1),'*')
%             hold off
   
   
   
    hold off
    h=figure(1);
        rect = [40, 1, 700, 270];
        frame=getframe(h,rect);
        imind=frame2im(frame);
        [imind,cm] = rgb2ind(imind,256);
    if simN==1
        imwrite(imind,cm,'1.gif','gif', 'Loopcount',inf,'DelayTime',0.01);%第一次必须创建！
    else
        imwrite(imind,cm,'1.gif','gif','WriteMode','append','DelayTime',0.01);
    end  
figure(11);
    width=800;%宽度，像素数
    height=270;%高度
    left=635;%距屏幕左下角水平距离
    bottem=200;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height],'color','w')
plot(track.outer(1,:),track.outer(2,:),'r','linewidth',3)
%     axis([-50,50,25,50,]);
    hold on
    plot(track.inner(1,:),track.inner(2,:),'r','linewidth',3)
%     plot(track2.outer(1,:),track2.outer(2,:),'k')
%     plot(track2.inner(1,:),track2.inner(2,:),'k')
%     plot(b(:,1),b(:,2),'g.')
%     plot(b(:,3),b(:,4),'g.')
    plot(ppval(traj.ppx,mod(x(7,:),tl)),ppval(traj.ppy,mod(x(7,:),tl)),':k')
    plot(x(1,:),x(2,:),'b')
    carBox(x(:,1),ModelParams.W/2,ModelParams.L/2)
    if ~isempty(Y)
        for i=1:size(Y,2)
            carBox_1(Y(:,i),ModelParams.W/2,ModelParams.L/2)
        end
    end
    xlabel('X [m]')
    ylabel('Y [m]')
    axis([0+x(7,1)*0,100*0+x(7,1)*0+73,-35,10]); 
%    axis([0+x(7,1)*0-3,100*0+x(7,1)*0+150,28,48,]);
%     axis([0+x(7,1)-10,100+x(7,1),25,50,]);
                                                       h=figure(11);
                                                        rect = [40, 1, 700, 270];
                                                        frame=getframe(h,rect);
                                                        imind=frame2im(frame);
                                                        [imind,cm] = rgb2ind(imind,256);
                                                    if simN==1
                                                        imwrite(imind,cm,'2.gif','gif', 'Loopcount',inf,'DelayTime',0.01);%第一次必须创建！
                                                    else
                                                        imwrite(imind,cm,'2.gif','gif','WriteMode','append','DelayTime',0.01);
                                                    end     
   
   
    figure(2)
    hold on
    width=600;%宽度，像素数
    height=500;%高度
    left=1300;%距屏幕左下角水平距离
    bottem=500;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height],'color','w')
    subplot(3,1,1)
    hold on
    plot(([0:N-1]+simN)*Ts,u(1,:))
    xlabel('time [s]')
    ylabel('D [-]')
    subplot(3,1,2)
    hold on
    plot(([0:N-1]+simN)*Ts,u(2,:))
    xlabel('time [s]')
    ylabel('\delta [rad]')
    subplot(3,1,3)
    hold on
    plot(([0:N-1]+simN)*Ts,u(3,:))
    xlabel('time [s]')
    ylabel('v_{\Theta} [m/s]')
figure(6666)
    
    width=550;%宽度，像素数
    height=200;%高度
    left=1300;%距屏幕左下角水平距离
    bottem=500;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height],'color','w')
        hold on
%     plot(([-1:50-2]+simN)*Ts,u(2,1:50)*180/pi)
    plot(([-1:N-2]+simN)*Ts,u(2,:)*180/pi)
    xlabel('time [s]')
    ylabel('\delta [deg]')
    
    figure(3)
    subplot(3,1,1)
    plot([0:N]*Ts,x(3,:))
    xlabel('time [s]')
    ylabel('\phi [rad]')
    subplot(3,1,2)
    plot([0:N]*Ts,x(6,:))
    xlabel('time [s]')
    ylabel('\omega [rad/s]')
    subplot(3,1,3)
    plot([0:N]*Ts,x(7,:))
    xlabel('time [s]')
    ylabel('\Theta [m]')
    
    figure(4)
     width=600;%宽度，像素数
        height=400;%高度
        left=1250;%距屏幕左下角水平距离
        bottem=0;%距屏幕左下角垂直距离
        set(gcf,'position',[left,bottem,width,height],'color','w')
    subplot(3,1,1)
    plot([0:N]*Ts,x(3,:))
    xlabel('time [s]')
    ylabel('YawAngle [rad]')
    subplot(3,1,2)
    plot([0:N]*Ts,x(4,:))
    xlabel('time [s]')
    ylabel('Vx [m/s]')
    subplot(3,1,3)
    plot([0:N]*Ts,x(5,:))
    xlabel('time [s]')
    ylabel('Vy [m/s]')

    pause(simTime)

end

