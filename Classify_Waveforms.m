function [class,Ptail,PP,RatioIP,kurt,PPloc] = Classify_Waveforms(SAT,WD,sigma0,classIN)

%CLASSIFY_WAVEFORMS adds the classes 'sea ice', 'specular', and 'undefined'
%to classIN (read from L1b file). The classification is applied to
%normalized waveforms and based on sigma0 and/or several waveform shape
%parameters

%The classification criteria are emperically determined (using among others
%optical imagery dat a), as well as on:
% -Poisson, J. C., Quartly, G. D., Kurekin, A. A., Thibaut, P., Hoang, D.,
%  & Nencioli, F. (2018). Development of an ENVISAT altimetry processor
%  providing sea level continuity between open ocean and Arctic leads.
%  IEEE Transactions on Geoscience and Remote Sensing, 56(9), 5299-5319.
% -Schulz, A. T., & Naeije, M. (2018). SAR Retracking in the Arctic:
%  Development of a year-round retrackers system. Advances in Space
%  Research.
% -Wernecke, A., & Kaleschke, L. (2015). Lead detection in Arctic sea ice
%  from CryoSat-2: quality assessment, lead area fraction and width
%  distribution. The Cryosphere, 9(5), 1955-1968.

%Input:
%SAT:     satellite mission from which data are processed ('CS'/'S3A'/'S3B')
%WD:      normalized waveform(s)
%sigma0:  sigma naught values in dB
%classIN: classification results read from L1b file

%Output:
%class:  0/1/2/3/4/5/99; where:
%   0 = open ocean or semi-enclosed seas; 
%   1 = enclosed seas or lakes; 
%   2 = continental ice; 
%   3 = land;
%   4 = specular; 
%   5 = sea ice;
%   6 = land-contaminated waveform; (determined by quality of fit)
%  99 = undefined.

%% Preliminaries
%Nr of bins/samples in any waveform & Nr of waveforms
[NrBins,NrWDs] = size(WD);

%Set width used to compute local pulse peakiness
switch SAT
    case {'CS','CryoSat'}
        width  = 3;
        tail   = 100;
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
        width  = 5;
        tail   = 35;
        sigma0 = sigma0 + 10;
    otherwise
        error('SAT: %s not recognized',SAT)
end

%% Compute waveform parameters
%Pulse Peakiness (PP) is the ratio of the maximum power and the accumulated
%signal power (first defined by Laxon (1994))
PP   = max(WD)./sum(WD)*(NrBins/256);

%Index of largest peak
[~,IDXpeak] = max(WD);

%Sample kurtosis of waveforms
kurt = kurtosis(WD);

%Determine local PP, Nr of peaks, start of leading edge, and power in tail
[PPloc,NrPeaks,StartLE,Ptail,RatioIP] = deal(nan(1,NrWDs));
for i = 1:NrWDs
    %Find start of leading edge
    StartLE(i) = find(WD(:,i) > 0.02,1);
    
    %Local PP is the peakiness of a part of the waveform; here, the part
    %surrounding the maximum peak. Local PP is used to distinguish between
    %true specular waveforms and specular waveforms with a tail
    if IDXpeak(i) > width && IDXpeak(i) < NrBins-width
        PPloc(i) = sum(WD(IDXpeak(i)-width:IDXpeak(i)+width,i))/sum(WD(:,i));
    else
        PPloc(i) = 10;
    end
    
    %Assess whether there are more significant peaks besides the main peak.
    %The significance of peaks is assessed based on their prominence and
    %distance to the main peak.
    NrPeaks(i) = numel(findpeaks(smooth(WD(:,i))./max(smooth(WD(:,i))),1:NrBins,'MinPeakProminence',0.05,'MinPeakDistance',5));
    
    %Capture additional peaks preceding the main peak.
    % if NrPeaks(i) == 1 & StartLE(i) < 150 & IDXpeak(i) <= 250, NrPeaks(i) = numel(findpeaks(WD(1:IDXpeak(i)+6,i),1:IDXpeak(i)+6,'MinPeakHeight',0.5,'MinPeakProminence',0.05)); end

    %Sum of total power in the tail (x no. of bins after the main peak).
    %The value has to be neglectable for lead waveforms and < 8 for ocean
    %ones.
    Ptail(i)   = nansum(WD(IDXpeak(i)+tail:end,i));
    
    %Ratio between the integrated power acquired before and after the max peak.
    %A ratio of ~0.5 indicates a specular waveform (leads). A ratio > 1 indicates
    %significant power before the max peak (sea ice/land contamination)
    RatioIP(i) = sum(WD(1:IDXpeak(i)-2,i))/sum(WD(IDXpeak(i)+2:end,i));
end

%% Apply classification
class = classIN*NaN;
%Ocean
class(PP < 0.05 & (NrPeaks == 1 & kurt < 15) & PPloc < 0.45 & ((sigma0 > 0 & sigma0 < 15) | isnan(sigma0)) & Ptail > 1 & Ptail < 8 & class ~= 0) = 0;  
%Lead/calm waters
class(StartLE > 30 & PP > 0.25 & NrPeaks == 1 & Ptail < 0.5 & (sigma0 > 25 | isnan(sigma0))) = 3;  %RatioIP < 0.65 & 
%Ice/land
class((NrPeaks > 1 | kurt > 12) & PP < 0.15 & class ~= 2 & class ~= 3) = 1;
%Set class to undefined for waveforms whose leading edge starts after bin 150
class(StartLE > 150*(NrBins/256) | RatioIP == Inf) = 99;

% [idx,C] = kmeans([PP;NrPeaks;kurt;PPloc;RatioIP;sigma0;Ptail;StartLE;IDXpeak]',4,'Display','iter','Replicates',5);
% [idx,C] = kmedoids([PP;double(NrPeaks > 1);kurt;PPloc;RatioIP;sigma0;Ptail;StartLE;IDXpeak]',4,'Replicates',5);
% idx = dbscan(X,epsilon,minpts);

end
