clear;
clc;
scenario='circle';   %'straight';% %   %%%%% 

switch scenario
    case 'circle'
%         cellPatient = xlsread('track.xlsx','inner','B1:NN6'); 
%         Trackdata.inner(1,:)=cellPatient(1,:)-0.6*0 ;
%         Trackdata.outer(1,:)=cellPatient(3,:)-0.6*0;
%         Trackdata.center(1,:)=cellPatient(5,:)-0.6*0 ;
%         Trackdata.inner(2,:)=cellPatient(2,:)+0.3*0 ;
%         Trackdata.outer(2,:)=cellPatient(4,:)+0.3*0 ;
%         Trackdata.center(2,:)=cellPatient(6,:)+0.3*0 ;
        cellPatient = xlsread('track.xlsx','inner','B29:FB34'); 
        Trackdata.inner(1,:)=cellPatient(1,:) ;
        Trackdata.outer(1,:)=cellPatient(3,:);
        Trackdata.center(1,:)=cellPatient(5,:) ;
        Trackdata.inner(2,:)=cellPatient(2,:)+0.3 ;
        Trackdata.outer(2,:)=cellPatient(4,:)+0.3 ;
        Trackdata.center(2,:)=cellPatient(6,:)+0.3 ;
    case 'straight'    
        cellPatient = xlsread('track.xlsx','B17:ABW22'); 
        Trackdata.inner=cellPatient(1:2,:) ;
        Trackdata.outer=cellPatient(3:4,:) ;
        Trackdata.center=cellPatient(5:6,:) ;
end

save Trackdata
clear
