function class = Classify_Waveform(unNORMdata, NORMfactor, sigma0, method)

%CLASSIFY_WAVEFORMS classifies normalized waveforms into ice, lead, ocean,
%or undefined waveforms based on sigma0 and/or several waveform parameters

%The classification criteria are emperically determined (using among others
%optical imagery data), as well as on:
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
%WD:     normalized waveform(s)
%sigma0: sigma naught values in dB

%Output:
%class:  0/1/2/3; where 1 = ice/land, 2 = ocean, 3 = lead,  0 = undefined

%% Preliminaries
NORMdata   = NORMfactor .* unNORMdata;
%Nr of bins/samples in any waveform
WD = NORMdata; 
NrBins = size(WD,1);  
MAX = max(unNORMdata) ; 

%Nr of waveforms
NrWDs  = size(WD,2);

%Set width, used to compute local pulse peakiness
if NrBins == 128
    width = 3;
elseif NrBins == 256
    width = 5;
end

%% Compute waveform parameters
%Pulse Peakiness (PP) is the ratio of the maximum power and the accumulated
%signal power (first defined by Laxon (1994))
PP     = max(WD)./sum(WD);

%Index of largest peak
[~,IDXpeak] = max(WD);

%Sample kurtosis of waveforms
kurt = kurtosis(WD);

%Determine local PP, Nr of peaks, start of leading edge, and power in tail
[PPloc,NrPeaks,StartLE,Ptail] = deal(nan(1,NrWDs));

for i = 1:NrWDs %NrWDs
    %Local PP is the peakiness of a part of the waveform; here, the part
    %surrounding the maximum peak. Local PP is used to distinguish between true
    %specular waveforms and specular waveforms with a tail
    if IDXpeak(i) > width && IDXpeak(i) < NrBins-width
        PPloc(i) = sum(WD(IDXpeak(i)-width:IDXpeak(i)+width,i))/sum(WD(:,i));
    else
        PPloc(i) = 10;
    end
    
    %Assess whether there are more significant peaks beside the main peak. The
    %significance of peaks is assessed based on their relative height and
    %prominence, and distance to the main peak. Peaks that are further away
    %have to be less high.
    NrPeaks(i) = numel(findpeaks(max(0,filter(-smooth_diff(8),1,WD(:,i))),1:NrBins,'MinPeakProminence',0.01,'MinPeakDistance',10));
    
    %Capture additional peaks preceding the main peak.
    if NrPeaks(i) == 1 && IDXpeak(i)<126
        NrPeaks(i) = numel(findpeaks(WD(1:IDXpeak(i)+2,i),1:IDXpeak(i)+2,'MinPeakHeight',0.5,'MinPeakProminence',0.05)); 
    end

    %Sum of total power in the tail (bins after the main peak+3). The value 
    %has to be neglectable for lead waveforms, but significant for ocean
    %ones. 
    Ptail(i)   = sum(WD(IDXpeak(i)+3:end,i));
    
    %Find start of leading edge
    StartLE(i) = find(WD(:,i) > 0.05,1);
end


%% Apply classification
class = zeros(1,NrWDs);

switch method 
    
    case 'Inger'
        %Ocean
        class(PP < 0.065 & NrPeaks == 1 & PPloc < 0.3 & (Ptail > 20 & Ptail < 40) & ((sigma0 > 0 & sigma0 < 10) | isnan(sigma0)) ) = 2; 
        %Lead/calm waters
        %class(PP > 0.25 & NrPeaks == 1 & Ptail < 0.5 & (sigma0 > 25 | isnan(sigma0))) = 2;  
        %Ice/land
        %class((NrPeaks > 1 | kurt > 12) & PP < 0.15 & class == 0) = 1;
        %Lead/calm waters
        class(PP > 0.25 & NrPeaks == 1 ) = 3;  
        %Ice/land
        class((NrPeaks > 1 | kurt > 12) & PP < 0.15 & class == 0) = 1;
        %Set class to undefined for waveforms whose leading edge starts after bin 150
        class(StartLE > 150) = 0;
        
        
    case 'Rose'
        %Lead/calm waters
        class(PP > 0.25 & NrPeaks == 1 ) = 3;  
        %Ice/land
        class((NrPeaks > 1 | kurt > 12) & PP < 0.15 & class == 0) = 1;
        %Set class to undefined for waveforms whose leading edge starts after bin 150
        class(StartLE > 150) = 0;
        
    case 'Peacock'
        %Lead/calm waters
        class(PP > 1.8 ) = 3;  
        %Ice/land
        class( PP < 1.8 ) = 1;
        %Set class to undefined for waveforms whose leading edge starts after bin 150
        class(StartLE > 150) = 0;
        
        
    case 'Kmed'
        class = fun_Kmed_KNN_classification(unNORMdata); 
    % AMBIGUOUS = 0
    % ICE = 1
    % LEAD = 3
    % OCEAN = 2
    % OUTLIERS = 4
        
end