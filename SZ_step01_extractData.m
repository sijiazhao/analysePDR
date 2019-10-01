addpath('./Functions');
clear;close all;
exptID = 'Expt10';
ifNormalisation = 1; % 2017-05-27 Move normalisation (convert to z-score) before baseline correction

%% Load Experimental setup parameters
anlmode = 'PDRatOnset';
timelock_code = 1; % ONSET - 20170525
tw_epoch = [-.5 8]; % time window for epoch [s]
tw_bc = [-.5 0]; % time window for epoch baseline correction [s]

%% Define and make directories
path_root = 'E:\Sijia\EOA\';
path_in = [path_root,'Data\' exptID '\']; % where is the data?
path_out = [path_root,'Results\',exptID,'\']; mkdir(path_out); % where is the analysed result?

%% Parameters for the analysis
ifrejectbadtrials = 1;
ifinterpolation = 0; %20161006

%% Retrieve subject informtion from excel file
[~,~,expMat] = xlsread([path_root,'Data\','subject_EYE_' exptID '.xlsx']);
sublist = expMat(:,1);  %Column 1: subject index
FsMatrix = expMat(:,4); %Column 5: sampling frequency, you can set it as 1000Hz
neyes = expMat(:,5);    %Column 6: how many eyes recorded? 1 or 2?
chooseeyes = cell2mat(expMat(:,6));     %Column 7: which eye to keep?

%% --

blocklist = 1:4;
% ntrial = 4*32;
disp(['# of blocks: ' numel(blocklist)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------- Start to convert -------------
% trialList = [];
P0 = {}; %pupil
summary.datalost = {};
summary.chooseeye = {};
mkdir([path_out 'raw']);

for subj = 1:length(sublist)
    disp(['**** Subject ' sublist{subj} ' conversion begins ****']);
    
    smpfreq = FsMatrix{subj};
    timeaxis = (tw_epoch(1)*smpfreq:tw_epoch(2)*smpfreq)/smpfreq;
    
    filename_raw = [path_out 'raw/raw_' sublist{subj}];
    if ~exist([filename_raw '.mat'],'file')
        
        for block = blocklist
            %% Load experimental parameters
            
            switch exptID
                case 'Expt10'
                    EyelinkName = [path_in, 'c' sublist{subj} '_b' num2str(block)]; %% {MySpecific} Adjust for your own experiments
                    filename_results = [path_in 'result_c' sublist{subj} '_b' num2str(block) '.mat'];
                otherwise
                    EyelinkName = [path_in, sublist{subj} '_b' num2str(block)]; %% {MySpecific} Adjust for your own experiments
                    filename_results = [path_in 'result_' sublist{subj} '_b' num2str(block) '.mat'];
            end
            load(filename_results, '-mat'); %e.g. stim_cS4_b6.mat
            
            %% Sc is a structure storing details about the screen.
            % The following three parameters are essential to compute visual degree.
            % In the present script, the visual degree is used to detect microssaccades.
            Sc.width =  ExpCond.ScWidth; % Screen width
            Sc.dist2Sbj = ExpCond.distSc_Sbj; % Disctance between screen and subject
            Sc.rect = ExpCond.rect; % Monitor resolution
            
            %% Load pupil data (* requires asc2data)
            
            [eyedata,pupildata,time,starttime,smpfreq] = asc2data([EyelinkName,'_sample.asc'],smpfreq,'binoc');
            raw.eyedata{block} = eyedata'; %LX,LY,RX,RY
            raw.pupildata{block} = pupildata'; %pupil dimeter of L and R
            raw.time{block} = time;
            raw.starttime{block} = starttime;
            raw.smpfreq{block} = smpfreq;
            
            %% Load Onsets of Trials
            filename = [EyelinkName, '_event.asc'];
            identifier = '%s %s %s %s'; % myspecific!  % (e.g.) MSG	1020407 Trial:2 ./StimData/S3X/2_25_CD_4_s1_250.wav
            fid = fopen(filename);
            tmp = textscan(fid, identifier);
            fclose(fid);
            
            ELindex = find(strcmp('MSG',tmp{1})==1 & strncmp('Trial',tmp{3},5)==1 & strncmp('Onset',tmp{4},5)==1); % tmp{4} stores the trial name
            
            if length(ELindex)==ExpCond.numTrial
                disp('Eyelink messages have been correctly detected')
                for k=1:ExpCond.numTrial
                    ELtime(block,k) = str2double(char(tmp{2}(ELindex(k))))-starttime;
                end
            else
                disp('! Some of Eyelink Messages are Missing !')
                
                tmptime = zeros(1,1);
                for k=1:length(ELindex)
                    tmptime(k) = str2double(char(tmp{2}(ELindex(k))))-starttime;
                end
                % CAN NOT ADOPT WHEN THE 1ST MESSAGE IS MISSING
                a=2;
                ELtime(block,1) = tmptime(1);
                for k=2:ExpCond.numTrial
                    if abs(tmptime(k)-tmptime(k-1)-mode(diff(tmptime)))<5
                        ELtime(block,k) = tmptime(k);
                    else
                        ELtime(block,k) = tmptime(k-1)+mode(diff(tmptime));
                        tmptime = [tmptime(1:k-1) ELtime(block,k) tmptime(k:end)];
                    end
                end
            end
            ELtime(block,:) = ELtime(block,:) - floor(ExpCond.audlag);
            
            
        end %block ends
        save([filename_raw '.mat'],'raw','ELtime');
    else
        load([filename_raw '.mat'],'raw','ELtime');
    end
    
    %% ---- %%
    disp(['**** Subject ' sublist{subj} ' : Epoch begins ****']);
    
    ps = [];
    %     bs = [];
    %     gs = [];
    %     mss = [];
    pds = [];
    pcs = [];
    pdds = [];
    pcds = [];
    summary.datalost{subj} =[];
    summary.chooseeye{subj} = [];
    for block = blocklist
        
        %% (prep.) Epoch pupildata and eyedata
        p = []; g = []; b = []; ms = []; pdEvent = []; pcEvent = [];
        eyedata     = raw.eyedata{block}; %LX,LY,RX,RY
        pupildata   = raw.pupildata{block}; %pupil dimeter of L and R
        time        = raw.time{block};
        starttime   = raw.starttime{block};
        smpfreq     = raw.smpfreq{block};
        eltime      = ELtime(block,:)'; %EyeLink 'time' points
        
        %% (prep.) MS
        mymsonset = MSraw{block}.left.onset;
        myms = zeros(1,length(time));
        myms(mymsonset) = 1;
        
        %% Handle blinks
        blink = zeros(length(time),1);
        blink(pupildata(:,chooseeyes(subj))==0) = 1;
        pupildata = ReplaceBlinksWithZeros(pupildata);
        for ip = 1:size(pupildata,2)
            Pn = pupildata(:,ip);
            eyedata(Pn==0,ip:ip+2) = NaN;
            Pn(Pn == 0) = NaN;
            [Pn,~]=foutlierremoval(Pn,2);
            pupildata(:,ip) = Pn;
        end
        
        %% Handle half blinks (based on pupil diameter derivative)
        tx = linspace(1,length(time)/smpfreq,length(time));
        eyename = {'Left','Right'};
        
        pupildata0 = pupildata;
        
        figure(subj);clf;
        for ei = 1:2
            
            g1 = eyedata(:,(ei-1)*2+1);
            g2 = eyedata(:,(ei-1)*2+2);
            
            %             d = [nan; diff(g)];
            data = g1;
            rg11 = nanmedian(data) - 2*nanstd(data);
            rg12 = nanmedian(data) + 3*nanstd(data);
            Ir = find(data<rg11 | data>rg12);
            fprintf('..... Data removed based on gaze : %.2f\n', numel(Ir)/length(tx));
            pupildata(Ir,ei) = NaN;
            eyedata(Ir,(ei-1)*2+1) = NaN;
            
            
            data = g2;
            rg21 = nanmedian(data) - 2*nanstd(data);
            rg22 = nanmedian(data) + 2*nanstd(data);
            Ir = find(data<rg21 | data>rg22);
            fprintf('..... Data removed based on gaze : %.2f\n', numel(Ir)/length(tx));
            pupildata(Ir,ei) = NaN;
            eyedata(Ir,(ei-1)*2+2) = NaN;
            
            subplot(3,2,ei);
            hold on;
            plot(tx,g1,'Color',[.8 .8 .8]);
            plot(tx,g2,'Color',[.8 .8 .8]);
            xl = xlim;
            
            line(xl,[rg11 rg11],'Color','r');
            line(xl,[rg12 rg12],'Color','r');
            line(xl,[rg21 rg21],'Color','r');
            line(xl,[rg22 rg22],'Color','r');
            
            g1 = eyedata(:,(ei-1)*2+1);
            g2 = eyedata(:,(ei-1)*2+2);
            plot(tx,g1,'Color','k');
            plot(tx,g2,'Color','k');
            
            ylabel('gaze');
            title(eyename{ei});
            
            pn = pupildata(:,ei);
            d = [nan;diff(pn)];
            data = d;
            if max(abs(data)) >= 20
                rg1 = nanmedian(data) - nanstd(data);
                rg2 = nanmedian(data) + nanstd(data);
                Ir = find(data<rg1 | data>rg2);
                fprintf('..... Data removed based on pupil diameter derivative : %.2f\n', numel(Ir)/length(tx));
                pupildata(Ir,ei) = NaN;
                eyedata(Ir,ei:ei+1) = NaN;
            else
                rg1 = nanmedian(data) - 3*nanstd(data);
                rg2 = nanmedian(data) + 3*nanstd(data);
                Ir = find(data<rg1 | data>rg2);
                fprintf('..... Data removed based on pupil diameter derivative : %.2f\n', numel(Ir)/length(tx));
                pupildata(Ir,ei) = NaN;
                eyedata(Ir,ei:ei+1) = NaN;
            end
            
            subplot(3,2,2+ei);
            pn2 = pupildata(:,ei);
            hold on;
            plot(tx,pn,'Color',[.8 .8 .8]);
            plot(tx,pn2,'--k');
            hold off;
            ylabel('Pupil diameter [a.u.]');
            
            subplot(3,2,4+ei);
            hold on;
            d2 = [nan;diff(pn2)];
            plot(tx,d,'Color',[.8 .8 .8]);
            plot(tx,d2,'--k');
            xl = xlim;
            line(xl,[rg1 rg1],'Color','r');
            line(xl,[rg2 rg2],'Color','r');
            hold off;
            ylabel('Pupil diameter derivative [a.u.]');
        end
        
        
        datalost = [];
        for ei = 1:2
            a1 = sum(isnan(pupildata(:,ei)))/length(time);
            a2 = sum(isnan(pupildata0(:,ei)))/length(time);
            datalost(ei) = a2;
        end
        
        %         [datalost,chooseeye] = min(datalost); % choose the eye which has least data lost
        chooseeye = 1;
        datalost = datalost(chooseeye);
        fprintf('** %s block%d choose %s eye, which had %.2f data loss.\n',sublist{subj},block,eyename{chooseeye},datalost);
        
        summary.datalost{subj} = [summary.datalost{subj},datalost];
        summary.chooseeye{subj} = [summary.chooseeye{subj},chooseeye];
        
        %% Epoch
        t = pupildata(:,1)/smpfreq;
        
        epoch_idx_sampling = ceil(tw_epoch(1)*smpfreq):ceil(tw_epoch(2)*smpfreq); epoch_idx_sampling = epoch_idx_sampling';
        epoch_time_sampling = epoch_idx_sampling/smpfreq; epoch_time_sampling = epoch_time_sampling';
        idx_bc = find(epoch_time_sampling >= tw_bc(1) & epoch_time_sampling <= tw_bc(2));
        
        for trial = 1:numel(eltime)
            thisepoch = epoch_idx_sampling + double(eltime(trial));
            
            p(trial,:) = pupildata(thisepoch,chooseeye);
            for i = 1:size(eyedata,2)
                g(trial,:,i) = eyedata(thisepoch,i);
            end
            b(trial,:) = blink(thisepoch,1);
            ms(trial,:) = myms(thisepoch);
            
            pdEvent(trial,:) = PDEvent(thisepoch,chooseeye);
            pdDur(trial,:) = PDDur(thisepoch,chooseeye);
            
            pcEvent(trial,:) = PCEvent(thisepoch,chooseeye);
            pcDur(trial,:) = PCDur(thisepoch,chooseeye);
        end
        
        %% Normalised z-score
        if ifNormalisation
            p = (p-nanmean(reshape(p,[numel(p),1])))/nanstd(reshape(p,[numel(p),1])); % Normalised to z score
        end
        
        %% Baseline correction
        for trial = 1:size(p,1) % trial by trial
            Pn = p(trial,:);
            Pn = Pn - nanmedian(Pn(idx_bc));
            p(trial,:) = Pn;
        end
        
        ps = [ps;p];
        
    end %block ends
    
    P0{subj} = ps;
end
writefile = [path_out,'pupil_',anlmode,'.mat'];
save(writefile,'sublist','P0','summary','timeaxis','-v7.3');
disp('--- Step 1: Complete ---');
