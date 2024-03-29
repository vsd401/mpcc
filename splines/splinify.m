

function [traj, borders] = splinify(track)
%[traj, borders] = splinify(track)
%calculate the required spline and derivatives

Tr=track.center;%(track.outer+track.inner)/2;
Xt=Tr(1,:);
Yt=Tr(2,:);
% Tr = track.center;
% Xt = track.center(1,:);
% Yt = track.center(2,:);
[traj.ppx, traj.ppy err]=normalizedSplineInterp(Xt,Yt,1,'y'); 
%traj.ppx = spline(1:length(Xt),Xt);

%calculate derivatives of desired trajectory
traj.dppx=getSplineDerivatives(traj.ppx);                                  %求一阶微分
traj.dppy=getSplineDerivatives(traj.ppy);
traj.ddppx=getSplineDerivatives(traj.dppx);                                %求二阶微分
traj.ddppy=getSplineDerivatives(traj.dppy);

%compute borders of track spline approximation
[borders.pplx,borders.pply]=computeCenter(traj.ppx,traj.ppy,track.inner(1,:),track.inner(2,:),1);%%路边界的样条系数
[borders.pprx,borders.ppry]=computeCenter(traj.ppx,traj.ppy,track.outer(1,:),track.outer(2,:),1);

% [borders.pplx,borders.pply]=computeCenter(traj.ppx,traj.ppy,track.inner(1,:),track.inner(2,:),1);
% [borders.pprx,borders.ppry]=computeCenter(traj.ppx,traj.ppy,track.outer(1,:),track.outer(2,:),1);


borders.dpplx=getSplineDerivatives(borders.pplx);
borders.dpply=getSplineDerivatives(borders.pply);
borders.dpprx=getSplineDerivatives(borders.pprx);
borders.dppry=getSplineDerivatives(borders.ppry);

%compute center (for compatibility and tests)
[borders.ppcx,borders.ppcy]=computeCenter(traj.ppx,traj.ppy,Tr(1,:),Tr(2,:),1);
borders.dppcx=getSplineDerivatives(borders.ppcx);
borders.dppcy=getSplineDerivatives(borders.ppcy);
borders.trackwidth=0.37;
end