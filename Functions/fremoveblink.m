function newPn = fremoveblink(Pn, smpfreq, remtime)
if nargin <3
    remtime = [-0.15 0.15]; %unit:secs
end

u  = [0; diff(Pn<=0)];
u1 = find(u==1);       % blink onsets
u2 = find(u==-1);      % blink offsets
if numel(u1) > numel(u2) && u1(end)>u2(end)
    u2 = [u2; numel(Pn)];
elseif numel(u1) < numel(u2) && u2(1)<u1(1)
    u1 = [1; u1];
end
v = [u1+remtime(1)*smpfreq, u2+remtime(2)*smpfreq]; % before to after blink
v(find(v<1)) = 1;
v(find(v>numel(Pn))) = numel(Pn);

newPn = Pn;
for i = 1:size(v,1)
    newPn(v(i,1):v(i,2))=0;
end
end