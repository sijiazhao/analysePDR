function [pt1,pt2] = findpeakshistogram(xxx,edges)
%special for michelle
if min(size(xxx)) ~= 1, error('matrix xxx should be a vector'); end
bins = [];
for k = 1:numel(edges)-1
    thisbin = xxx(find(xxx >= edges(k) & xxx < edges(k+1)));
    bins = [bins, numel(thisbin)];
end
[maxv,maxt] = max(bins);
% maxt = mean(edges([maxt,maxt+1]));
pt1 = edges(maxt); %peak bin start
pt2 = edges(maxt+1); %peak bin end
end
