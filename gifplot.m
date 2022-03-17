clear
x0=5*0;
y0=5*0;
a=8;
b=6;
t=linspace(0,2*pi,100);
AA=a*sin(t)-x0;
BB=b*cos(t)-y0;
figure(1)
plot(AA,BB)


close all

% % % % % % % % t=linspace(0,2*pi,100);
% % % % % % % % x=10*sin(t)-5;
% % % % % % % % y=5*cos(t)-5;
% % % % % % % % phi=pi/4;
% % % % % % % % A=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[x;y];
% % % % % % % % plot(A(1,:),A(2,:))
% % % % % % % % hold on
% % % % % % % % 
% % % % % % % % x=10*sin(t)*0.5-5;
% % % % % % % % y=5*cos(t)*0.5-5;
% % % % % % % % A=[cos(phi) -sin(phi);sin(phi) cos(phi)]*[x;y];
% % % % % % % % plot(A(1,:),A(2,:))




% t=0:0.02:10;  Nt=size(t,2);
% x = cos(2*t).*(cos(t).^2);
% y = sin(2*t).*(sin(t).^2);
% for i=1:Nt;
% cla;hold on;
% figure(11)
% plot(x,y)
% plot(x(i),y(i),'o');
% h=figure(11);
% frame=getframe(h);
% imind=frame2im(frame);
% [imind,cm] = rgb2ind(imind,256);
%     if i==1
%         imwrite(imind,cm,'2.gif','gif', 'Loopcount',inf,'DelayTime',0.01);%第一次必须创建！
%     else
%         imwrite(imind,cm,'2.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
% end

% for i = 1:simN
%     figure(8)
%     plot(i,i^2,'*')
% 
% 	h = figure(8);
% 	% 下面一行替换为想要显示的内容
% 
% 	imshow([num2str(i),'.jpg']);
% 
% 
%     frame = getframe(h);
% 
%     im = frame2im(frame);
% 
%     [imind,cm] = rgb2ind(im,256);
% 
%  if i == 1
% 	imwrite(imind,cm,[num2str(i),'.gif'],'gif', 'Loopcount',inf,'DelayTime',0.01);
%  else
% 	imwrite(imind,cm,[num2str(i),'.gif'],'gif','WriteMode','append','DelayTime',0.01);
%  end
% end