

function l=splinelength(ppx, ppy,umin,umax)
%l=splinelength(ppx, ppy,umin,umax)
%compute the length of a 2D spline between parameter umin and umax
%ppx and ppy are spline parameter. See spline function for more details
%umin and umax are vectors of the same length
%and for all i, umin(i)<=umax(i)
%checked
assert(length(umin)==length(umax));

dppx=getSplineDerivatives(ppx);
dppy=getSplineDerivatives(ppy);

h = @(u) sqrt(ppval(dppx,u).^2+ppval(dppy,u).^2);
l=zeros(length(umin),1);
for i=1:length(umin)
    l(i)=quad(h,umin(i),umax(i));                                         %积分函数，相当于在每一步长中把X和Y方向沿着切线方向积分得到每一个片段的长度，i.e.  compute the length of a 2D spline
end

end
