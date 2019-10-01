function [x,y]=fsubplotnum(n)
if n<=4 
    x=n; y=1;
else
    x=ceil(n/4);
    y=4;
end
    
end