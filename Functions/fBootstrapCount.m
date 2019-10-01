function [sig, pos, neg]=fBootstrapCount(x1,P,M)
% This function uses bootstrap statistics difference vector x1 and
% outputs time segments where two conditions are significantly different at
% an alpha value of P, and for at least M consecutive samples.
%
% INPUT:
%   - x1: matrix of bootstrap difference vector (TIME*BOOTSTRAP ITERATIONS)
%   - P: p value [default: 0.001]
%   - M: minimum number of consecutive significant points to define a
%        significant segment
%
% OUTPUT:
%   - sig: vector of same size as x1 with:
%             0 for non significant segments
%            +1 for positive significant difference
%            -1 negative significant difference
%   - pos: segments where difference is positive
%   - neg: segments where difference is negative
%
% Nicolas Barascud (Last edited 01/07/2015)

if nargin < 3; M = 20; end
if nargin < 2; P=0.001; end

[Nsamples, B]=size(x1);
clear pos neg
%% Take portions of signal w/ 95% of bootstrap resampling above or below 0
CI = [P 1-P]*B; % confidence interval

pos_num = (x1==abs(x1));
pos_num = sum(pos_num,2); % number of positive resampling

neg_num = pos_num < CI(1);
pos_num = pos_num > CI(2);

% Case where first or last number is 1, append 0
pos_num = [0; pos_num(:); 0];
neg_num = [0; neg_num(:); 0];

% Positive differences
diff_vect=diff(pos_num);
if sum(abs(diff_vect))~=0
    pos(:,1) = find(diff_vect > 0);
    pos(:,2) = find(diff_vect < 0)-1;
    pos(:,3) = pos(:,2)-pos(:,1)+1;
else
    pos = zeros(1,3);
end

% Negative differences
diff_vect=diff(neg_num);
if sum(abs(diff_vect))~=0
    neg(:,1) = find(diff_vect > 0);
    neg(:,2) = find(diff_vect < 0)-1;
    neg(:,3) = neg(:,2)-neg(:,1)+1;
else neg = zeros(1,3);
end

%% IMPORTANT : once resampling has been done to find new baseline, only
% take significant intervals of duration greater than chance

% Then, find the strings of zeroes with a duration greater than or equal to
% some value (such as 3):
% stringIndex = (duration >= 3);
% startIndex = startIndex(stringIndex);
% endIndex = endIndex(stringIndex);
segIdxPos = (pos(:,3) >= M);
segIdxNeg = (neg(:,3) >= M);

pos = pos(segIdxPos,:);
neg = neg(segIdxNeg,:);

% Build vector of zeroes, ones and minus ones of the same size as x1,
% indicating non-significant, positive and negative differences,
% respectively
sig = zeros(1,Nsamples);
try
    for i=1:size(pos,1)
        sig(pos(i,1):pos(i,2)) = 1;
    end
end
try
    for i=1:size(neg,1)
        sig(neg(i,1):neg(i,2)) = -1;
    end
end
%disp([num2str(LH(idx1)) 'ms (LH), ' num2str(RH(idx2)) 'ms (RH)'])

end
