function [pmEvent,pmAmpl,pmDur] = pmdetect(Pn,mode)
%% detect pupil dilation/constriction event based on the duration threshold
% Pn = pupil data (one column)
% mode = 'PD' or 'PC'
debugmode = 0;
[pks1,locs1,~,~] = findpeaks(Pn); %Peak
[pks2,locs2,~,~] = findpeaks(-Pn); %Trough
pks2 = -pks2;

% figure(1000);clf;
% plot(Pn);

% try
% trough followed by a peak

if numel(locs2)<numel(locs1)
    locs2 = [locs2, length(Pn)];
elseif numel(locs1)<numel(locs2)
    locs1 = [1, locs1];
end
locs = [locs1;locs2];

locsSign = [ones(size(locs1)); -1*ones(size(locs2))];
[locs,I] = sort(locs);
locsSign = locsSign(I);

% Remove troughs/peaks which are continuous
cc = find(diff(locsSign)==0); locsSign(cc)=[];locs(cc)=[];
if locsSign(1)>0; locs=[1, locs]; locsSign=[-1,locsSign]; end
if mod(numel(locs),2)~=0 % if number is not even
    if locsSign(1)<0; locs=[locs, numel(Pn)]; locsSign=[locsSign, 1]; end
end

locs2 = locs(find(locsSign == -1));
locs1 = locs(find(locsSign == 1)); %locs1 is peak, locs2 is trough


% Step 2: Find differential (concerning rates of change and slopes of curve)
x = 1:numel(Pn);
y = Pn;
dx = mean(diff(x)); % Find Mean Differece In �x� Values
dy = gradient(y,dx);

dypos = find(dy>0);
dyneg = find(dy<0);

% Find the time point of switch: shown as the points when green changes to red
dyy = zeros(size(dy));
dyy(dypos) = 1; switchToPos = find(diff(dyy) == 1);
dyy(dyneg) = -1; switchToNeg = find(diff(dyy) == -2);

if debugmode
    % ------
    figure(9999);clf; subplot(2,1,1); plot(Pn,'b');
    hold on; scatter(locs,Pn(locs),'r','d','filled', 'LineWidth',1.5);
    title('pupil size: peaks and troughs');
    subplot(2,1,2);plot(dy,'r');
    hold on;  scatter(dypos,dy(dypos),'r','o','filled', 'LineWidth',1.5);
    hold on;  scatter(dyneg,dy(dyneg),'g','o','filled', 'LineWidth',1.5);
    hold on;  scatter(switchToPos,dy(switchToPos),'b','d','filled', 'LineWidth',1.5);
    title('deviation of pupil change: red for dilation, green for constriction');
    % ------
end

%% Define 'approved' PD: threshold parameter for the pupil dilation length

if switchToNeg(1) < switchToPos(1)
    switchToPos = [1, switchToPos];
end
if numel(switchToNeg) < numel(switchToPos)
    switchToNeg = [switchToNeg, numel(Pn)];
end
%     switchToPos1 = [];
%     switchToNeg1 = [];
%     for i = 1:numel(switchToPos)
%         if (switchToNeg(i)-switchToPos(i))> thrs % the time between switch to neg (starts to constriction) and switch to pos (starts to dilate) needs to be longer than thrs
%             switchToNeg1 = [switchToNeg1; switchToNeg(i)];
%             switchToPos1 = [switchToPos1; switchToPos(i)];
%         end
%     end
%     switchToPos = switchToPos1;
%     switchToNeg = switchToNeg1;

%% --------------------------------------------------------------------
%    Check if switchToPos (switch point) falls into trough-peak pair
% and find the peak time for each onset time in order to compute the
% amplitude

switch mode
    case 'PD'
        la = switchToPos;
        lb = locs1;
        l1 = [];
        l2 = [];
        if ~isempty(la)
            for i = 1:numel(la)
                % find its peak time
                tmptmp = lb(find(lb > la(i))); %locs1 is peak, locs2 is trough
                if ~isempty(tmptmp) % there is a peak point after pd onset!
                    l2 = [l2; la(i)]; %pd onset
                    l1 = [l1; tmptmp(find(min(tmptmp-la(i))))]; %find the closest pd peak time
                end
            end
        end
        
        onsets = l2; % alsmots equal to switch to pos
        peaktimes = l1;
        % --------------------------------------------------------------------
        
        pmEvent = zeros(size(Pn));
        pmEvent(onsets) = 1;
        
        pmAmpl = NaN(size(Pn));
        pmAmpl(onsets) = Pn(peaktimes)-Pn(onsets);
        
        pmDur = NaN(size(Pn));
        pmDur(onsets) = peaktimes - onsets;
        
    case 'PC'
        
        la = switchToNeg;
        lb = locs2;
        l1 = [];
        l2 = [];
        if ~isempty(la)
            for i = 1:numel(la)
                % find its peak time
                tmptmp = lb(find(lb < la(i))); %locs1 is peak, locs2 is trough
                if ~isempty(tmptmp) % there is a peak point after pd onset!
                    l2 = [l2; la(i)]; %pd onset
                    l1 = [l1; tmptmp(find(min(la(i)-tmptmp)))]; %find the closest pd peak time
                end
            end
        end
        
        onsets = l2; % alsmots equal to switch to pos
        peaktimes = l1;
        % --------------------------------------------------------------------
        
        pmEvent = zeros(size(Pn));
        pmEvent(onsets) = 1;
        
        pmAmpl = NaN(size(Pn));
        pmAmpl(onsets)=Pn(onsets)-Pn(peaktimes);
        
        pmDur = NaN(size(Pn));
        pmDur(onsets) = peaktimes - onsets;
end
if debugmode
    %---------------
    figure(9999);clf; subplot(2,1,1); plot(Pn,'b');
    % hold on; scatter(locs1,pks1,'r','d','filled', 'LineWidth',1.5);
    % hold on; scatter(locs2,pks2,'g','d','filled', 'LineWidth',1.5);
    
    % hold on; scatter(locs,Pn(locs),'r','d','filled', 'LineWidth',1.5);
    
    hold on; scatter(onsets,Pn(onsets),'b','d','filled', 'LineWidth',1.5);
    hold on; scatter(peaktimes,Pn(peaktimes),'g','d','filled', 'LineWidth',1.5);
    
    subplot(2,1,2);plot(dy,'r');
    % hold on;  scatter(dypos,dy(dypos),'r','d','filled', 'LineWidth',1.5);
    % hold on;  scatter(dyneg,dy(dyneg),'g','d','filled', 'LineWidth',1.5);
    % hold on;  scatter(switchToPos,dy(switchToPos),'b','d','filled', 'LineWidth',1.5);
    hold on; scatter(onsets,dy(onsets),'b','d','filled', 'LineWidth',1.5);
    %---------------
end
% else
% catch
%     pmEvent = zeros(size(Pn));
%     pmAmpl = NaN(size(Pn));
%     pmDur = NaN(size(Pn));
% end
end
