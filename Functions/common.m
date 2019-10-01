function x = common(A,B)

tmpA = sort(A);
tmpA(diff(tmpA)==0)=[];
tmpB = sort(B);
tmpB(diff(tmpB)==0)=[];

tmp = sort([tmpA tmpB]);
x = tmp(diff(tmp)==0);