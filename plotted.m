figure(6667)
    width=550;%��ȣ�������
    height=200;%�߶�
    left=1300;%����Ļ���½�ˮƽ����
    bottem=500;%����Ļ���½Ǵ�ֱ����
    set(gcf,'position',[left,bottem,width,height],'color','w')
    plot([0:simN-1]*Ts,U_log(2,:)*180/pi,'r','linewidth',2)
    xlabel('time [s]')
    ylabel('\delta [deg]')
%     axis([0,5,-5 5])
    
figure(6668)
     width=550;%��ȣ�������
    height=200;%�߶�
    left=1300;%����Ļ���½�ˮƽ����
    bottem=500;%����Ļ���½Ǵ�ֱ����
    set(gcf,'position',[left,bottem,width,height],'color','w')
    plot([0:simN-1]*Ts,U_log(3,:),'r','linewidth',2)
    xlabel('time [s]')
    ylabel('v_x [m/s]')
    
    figure(6666)
    axis([0 8 -40 10])
    