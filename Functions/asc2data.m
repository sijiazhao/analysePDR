function [eyedata,pupil,time,starttime,smpfreq] = asc2data(filename,smpfreq,rectype,varargin)

%
%
%
%

if nargin==4
    method = varargin{1};
    if isnumeric(method)==1
        if length(method)==2
            medPos = method;
        else
            error('RECT should be 1x2 vector')
        end
    elseif ischar(method)==1
        if strcmp(method,'median')==1
            medPos = NaN;
            % not recommended
        else
            error('UNEXPECTED INPUT FOR RECT')
        end
    end
else
    medPos = [0 0];
end
dt = 1000/smpfreq;

tic;
identifier = '%u %s %s %f %s %s %s %s %s';
disp('LOADING FILE.....')
fid = fopen(filename);
tmp = textscan(fid, identifier);
fclose(fid);
disp([' - File load completed: ',num2str(toc),'[s]'])
time = tmp{1};
if isempty(find(diff(time)~=dt, 1))==0
%     startindex = 1 + find(diff(time)~=dt, 1, 'last' );
%     disp([' - DETECTED INVALID RECORDING!',char(10),' --- Starting sampling from ',num2str(startindex)]);
    dtime = diff(time);
    sx = 1 + find(dtime~=dt);
    jumptime = dtime(sx-1);
    disp(' - DETECTED INVALID RECORDING!');
    disp(['   total number of recording slip:',num2str(length(sx))]);
    for k=1:length(sx)
        disp(['[',num2str(k),'] sample:',num2str(sx(k)),' jump:',num2str(jumptime(k))])
    end
%     nt1 = input('STARTINDEX? (NONE: start from 1st sample)=');
%     nt2 = input('ENDINDEX? (NONE: end at last sample)=');
    nt1 = 1;
    nt2 = [];
    if isempty(nt1)==1
        startindex=1;
    else
        startindex = sx(nt1);
    end
    if isempty(nt2)==1
        endindex=length(time);
    else
        endindex = sx(nt2);
    end
else
    startindex = 1;
    endindex = length(time);
end
Ln = length(time);
time = time(startindex:endindex);
L = length(time);
starttime = time(1);

if strcmp(rectype,'monoc')==1
    pupil(1,1:L) = tmp{4}(startindex:endindex)';
    tmpdata = -30000*ones(2,Ln);
    nomissindex = cell(1,2);
    recnum=[2 3];
    p = 1;
elseif strcmp(rectype,'binoc')==1
    pupil = zeros(2,L);
    pupil(1,1:L) = tmp{4}(startindex:endindex);
    try
        pupil(2,1:L) = tmp{7}(startindex:endindex);
    catch err
       pupil(2,1:L) = sscanf(reshape([char(tmp{7}(startindex:endindex)) repmat(' ',length(tmp{7}(startindex:endindex)),1)]',1,[]),'%f')'; 
    end
    tmpdata =-30000*ones(4,Ln);
%     tmpdata =NaN*ones(4,Ln); %Sijia 2016/12/1
    nomissindex = cell(1,4);
    recnum=[2 3 5 6];
    p = 2;
end

a=0;
for i=recnum
    a = a+1;
    nomissindex{a} = find(strcmp('.',tmp{i})==0);
    tmpdata(a,nomissindex{a}) = sscanf(reshape([char(tmp{i}(nomissindex{a})) repmat(' ',length(nomissindex{a}),1)]',1,[]),'%f');
    if isnan(medPos)==1
        tmpdata(a,nomissindex{a}) = tmpdata(a,nomissindex{a})-median(tmpdata(a,nomissindex{a}));
    else
        tmpdata(a,nomissindex{a}) = tmpdata(a,nomissindex{a})-medPos(2-mod(a,2));
    end
end
disp('***** ALL PROCESS HAS BEEN COMPLETED *****')
toc
eyedata = tmpdata(1:a,startindex:endindex);
clear tmp tmpdata

% UPSAMPLING
if smpfreq==500
    time = bisample(time);
    for k=1:p
        tmp(k,1:L*2-1) = bisample(pupil(k,:));
    end
    pupil = tmp;
    tmp = zeros(a,L*2-1);
    for k=1:a
        tmp(k,:) = bisample(eyedata(k,:));
    end
    eyedata = tmp;
    smpfreq = 1000;
end