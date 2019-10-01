function s_boot = removeShortInterval(s_boot,M)
% remove the cluster whose length <M
% [~,tzero] = fFindClosest(timeaxis,0);
s_boot(find(s_boot==0)) = NaN;
s_boot(find(~isnan(s_boot))) = 1;
s_boot(find(isnan(s_boot))) = 0;

s = s_boot(1:length(s_boot));
sd = diff(s);
sd_start = find(sd==1); % start to be significant
sd_end = find(sd==-1); %end significance
if sd_end(1)<sd_start(1)
    sd_start = [1,sd_start];
end
if numel(sd_end) > numel(sd_start)
    sd_start = [1,sd_start];
end
if numel(sd_end) == numel(sd_start)-1
    sd_end = [sd_end,length(s_boot)];
end

sd_length = abs(sd_end-sd_start);

check = find(sd_length<M);
if ~isempty(check)
    for c = check
        the_naughty_cluster = (sd_start(c):1:sd_end(c));
        s_boot(the_naughty_cluster) = NaN;
    end
    sd_start(check) = [];
    sd_end(check) = [];
end

s_boot(find(s_boot==0)) = NaN;
