function pupil_out = ReplaceBlinksWithZeros(pupil_in,RemoveDur,CheckFlag)
% ReplaceBlinksWithZeros - Replace blinks, including transient portions,
% with zeros.
%Usage: pupil_out = ReplaceBlinksWithZeros(pupil_in,RemoveDur, CheckFlag)
%pupil_in: pupil matrix; Colmun1: time in ms; Colmun2-3: pupil data
%RemoveDur: Durations to be removed before and after the
%blinks in sec. (2-element vec; default [0.15, 0.15])
%CheckFlag: If non-zero, plot waveforms before and after the blink removal
%(default: 0)
%By SF, 2012/12/17

if nargin <2
    RemoveDur=[0.15, 0.15];
    if nargin <3
        CheckFlag=0;
    end
end
RemoveDurStart=RemoveDur(1); %Duration for removal before blink onset in sec.
RemoveDurEnd=RemoveDur(2); %Duration for removal after blink offset in sec.

%Output var
pupil_out=pupil_in;

%Extract data
t=pupil_in(:,1);
t=t-t(1);
NData=size(pupil_in,2)-1;

%sampling rate
myt=diff(t);
if myt(1) == 0, myt = find(myt~=0); end
Fs=1/myt(1)*1000;

%Indices for the analyses and removal
NRemoveStart=ceil(Fs*RemoveDurStart);
NRemoveEnd=ceil(Fs*RemoveDurEnd);

for iData=1:NData
    p=pupil_in(:,iData+1);
    
    %Detect blinks
    myp1=diff([p; p(end)]);
    myp2=diff([p(1); p]);
    IdxBlinkStart=find(myp2<0 & p==0); %Negavie slope -> zero
    IdxBlinkEnd=find(myp1>0 & p==0); %zero -> Positive slope
        
    if p(1) == 0, IdxBlinkStart = [1; IdxBlinkStart]; end
    
    %Replace the data with zeros for periods to be removed
    IdxRemoveStart=repmat(IdxBlinkStart,[1, NRemoveStart]) + repmat(-(1:NRemoveStart),[length(IdxBlinkStart),1]);
    IdxRemoveStart(IdxRemoveStart<1)=1;
    IdxRemoveEnd=repmat(IdxBlinkEnd,[1, NRemoveEnd]) + repmat(1:NRemoveEnd,[length(IdxBlinkEnd),1]);
    IdxRemoveEnd(IdxRemoveEnd>length(p))=length(p);
    
    if CheckFlag
        %Plot entire original waveform
        subplot(NData,1,iData)
        plot(t,p,'-')
    end
    
    %Replace
    p(IdxRemoveStart(:))=0;
    p(IdxRemoveEnd(:))=0;    
    
    pupil_out(:,iData+1)=p;
    
    if CheckFlag
        %Show the data after blinks removed.
        IdxUse=find(p>0);
        hold on
        plot(t(IdxUse),p(IdxUse),'r-')
        hold off
        
        xlabel('Time (sec)');
        ylabel('pupil diameter');
        title('Blue: Original; Red: After removal of blinks');
    end

end
