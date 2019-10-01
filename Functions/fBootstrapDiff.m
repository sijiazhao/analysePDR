function [mn,sd,all]=fBootstrapDiff(x1,x2,B,balanced,type)
% BOOTSTRAP resampling of difference vector between two populations.
%
% INPUTS
%	- x1 : matrix of observation for condition1 (TIME*REPETITIONS)
%	- x2 : matrix of observation for condition1 (TIME*REPETITIONS)
%	- B : number of bootstrap iterations
% - type : type of bootstrap to compute (see below)
%
% OUTPUTS
%	- mm : mean of bootstrap resampling
%	- sd : standard deviation of bootstrap resampling
% - all : matrix of all bootstrap iterations
%
% Last edited by Nicolas Barascud on 11/10/2011

if nargin < 5; type = 'between'; end
if nargin < 4; balanced = 1; end
if nargin < 3; B = 1000; end
if nargin < 2; error('Error: not enough arguments!'); end

if ~ismatrix(x1) || ~ismatrix(x2); error('Data must be at most 2D'); end

[Nsamples, Nrep] = size(x1);

% Creates bootstrap permutations in matrix p
% p will have dimensions B*Nsamples:
% --------------
% | 1 ... Nrep |
% | 1 ... Nrep |
% |    .       |- B
% |    .       |
% |    .       |
% | 1 ... Nrep |
% --------------
if balanced
    % balanced implementation ensures that all trial/repetitions have been
    % used at the end of the resampling
    p  = repmat(1:Nrep,B,1);
    p  = p(reshape(randperm(B*Nrep),B,Nrep)); % balanced implementation
    p2 = repmat(1:Nrep,B,1);
    p2 = p(reshape(randperm(B*Nrep),B,Nrep)); % balanced implementation
else
    p  = randi([1 Nrep],B,Nrep);
    p2 = randi([1 Nrep],B,Nrep);
end

all = zeros(Nsamples,B);
x = 0; h = waitbar(x,'Calculating Bootstrap');

if strcmp(type,'trialbased')
    % For each bootstrap iteration, uses different resamplings for both
    % conditions. For each bootstrap iteration, averages of conditions A
    % and B are computed from different (random) subject lists.
    % This is less robust to (e.g.) subject-specific DC shifts.
    for b = 1:B
        all(:,b) = mean(x1(:,p(b,:)),2) - mean(x2(:,p2(b,:)),2);
        x = b/B; waitbar(x); % update waiting bar
    end
    mn = mean(x1,2)-mean(x2,2);

elseif strcmp(type,'between')
    % For each bootstrap iteration, use the same random sample for both
    % conditions. For each bootstrap iteration, averages of conditions A
    % and B are computed from the SAME (random) subject list.
    % This is more robust to between-subject differences.
    for b = 1:B
        all(:,b) = mean(x1(:,p(b,:)),2) - mean(x2(:,p(b,:)),2);
        x = b/B; waitbar(x); % update waiting bar
    end
    mn = mean(x1,2) - mean(x2,2);

elseif strcmp(type,'diff')
    % This implementation uses the difference wave as input.
    % Not sure if this is different from second implementation
    x3 = x1-x2;
    for b = 1:B
        all(:,b) = mean(x3(:,p(b,:)),2);
        x = b/B; waitbar(x); % update waiting bar
    end
    mn = mean(x3,2);
end
close(h);

sd = std(all,1,2);
