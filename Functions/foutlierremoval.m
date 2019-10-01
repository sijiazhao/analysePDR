function [y,removedSubj]=foutlierremoval(x,numstd)
% function [y]=foutlierremoval(x,numstd)
% Updated on 22/12/2016
if nargin<2; numstd=2; end
meanRT = nanmean(x);
stdRT = nanstd(x);
try
    removedSubj=[find(x>(meanRT+numstd*stdRT)); find(x<(meanRT-numstd*stdRT))];
catch
    removedSubj=[find(x>(meanRT+numstd*stdRT)) find(x<(meanRT-numstd*stdRT))];
end
x(x>meanRT+numstd*stdRT) = NaN; x(x<meanRT-numstd*stdRT) = NaN;
y = x;
end