function y=makegrid(k,n,min,max)

if k==1
    y=linspace(min,max,n);  % equidistant grid
elseif k==2
    y=linspace(0,log(1+max-min),n); % exponential grid, more points at low assets, makes sense if curvature (or action) is there
    y=exp(y)-1;
    y=y+min;
else
    y=linspace(0,log(log(1+max-min)+1),n);  % exponential grid, even more points at low assets, same reason
    y=exp(exp(y)-1)-1;
    y=y+min;
end
end