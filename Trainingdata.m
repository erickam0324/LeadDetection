% Training data processing 
% Aquire all the data that has been labeled through the validation scheme 

clear all; 
%
D = '/Users/ericka/Desktop/Thesis/ValidationImageOutput';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:16%numel(N)
    LeadConfidence{ii} = load(fullfile(D,N{ii},'lead_confidence.mat')); % improve by specifying the file extension. 
    wfSA{ii} = load(fullfile(D,N{ii},'wfSA.mat'));
    sigma0{ii} = load(fullfile(D,N{ii},'sigma0.mat'));
end

% put all data together 
wf_all = []; 
conf_all = [];
sigma0_all = [];
for ii = 1:16%numel(N)
    wf_all = [wf_all wfSA{1,ii}.wftosave]; 
    sigma0_all = [sigma0_all sigma0{1,ii}.sigma0.']; 
    conf_all = [conf_all LeadConfidence{1,ii}.lead_confidence] ; 
end 

% Find the validation accuracy (require manual inspection)
% % input number of false ice and false leads seen by human eye : 
% FalseIce_val = [1,0,0,2,0,2,0] ; 
% FalseLead_val = [10, 1, 1, 13,4,4,8 ] ; 
% 
% for ii = 1:numel(N)
%     LeadConfidence{1,ii}.lead_confidence(LeadConfidence{1,ii}.lead_confidence ==4 ) =3 ; 
%     LeadConfidence{1,ii}.lead_confidence(isnan(LeadConfidence{1,ii}.lead_confidence)) = [] ; 
%     LeadConfidence{1,ii}.lead_confidence(LeadConfidence{1,ii}.lead_confidence == 2 ) =[] ; 
%     n_SA(ii) = numel(LeadConfidence{1,ii}.lead_confidence) ; 
%     Accuracy_val(ii) = 1 - (FalseIce_val(ii) + FalseLead_val(ii)) / (n_SA(ii)) ; 
%     TrueLead_val(ii) = sum(LeadConfidence{1,ii}.lead_confidence == 3) - FalseLead_val(ii); 
%     TrueIce_val (ii) = sum(LeadConfidence{1,ii}.lead_confidence == 1) - FalseIce_val(ii); 
%     TLR_val(ii) = TrueLead_val(ii) / (TrueLead_val(ii) + FalseIce_val(ii)) ; 
%     FLR_val(ii) = FalseLead_val(ii) / (FalseLead_val(ii) + TrueIce_val(ii)) ; 
% end
% 
% %find weighted average of the TLR and FLR of all the study areas 
% TLR_val_all = sum(TLR_val.* n_SA) / (sum(n_SA)) ; 
% FLR_val_all = sum(FLR_val.* n_SA) / (sum(n_SA)) ; 
% Accuract_val_all = sum(Accuracy_val.* n_SA) / (sum(n_SA)) ; 


%
% 
% wf_all = wfSA{1,4}.wftosave ; 
% conf_all =  LeadConfidence{1,4}.lead_confidence; 

% consider confidence '4' as leads too 
conf_all(conf_all ==4) = 3;

%normalize the data 
%Wf_Norm_Aw = 1;
%NORMfactor = 1./max(movmean(wf_all,Wf_Norm_Aw,1),[],1);
%wf_all   = NORMfactor .* wf_all;
normwf = wf_all - min(wf_all(:)) ; 
normwf = normwf ./ max(normwf(:)) ; 

% lets ignore ambiguous class for now 
wf_all(:,conf_all == 2)= [] ; 
normwf(:,conf_all == 2)= [] ; 
sigma0_all(conf_all == 2)= [] ; 
conf_all(conf_all==2) = [];

wf_all(:,isnan(conf_all))= [] ; 
normwf(:,isnan(conf_all))= [] ; 
sigma0_all(isnan(conf_all))= [] ; 
conf_all(isnan(conf_all)) = [];

nwf_all = numel(conf_all) ; 


% compute waveform parameters for each waveform 
[MAX_all, MAX_bin] = max(normwf) ; 
kurt_all = kurtosis(normwf) ; 
skew_all = skewness(normwf) ; 
PP_all = max(normwf)./sum(normwf);

% parameters from Inger : 

NrBins = size(normwf,1);  
%Set width, used to compute local pulse peakiness
if NrBins == 128
    width = 3;
elseif NrBins == 256
    width = 5;
end

[~,IDXpeak] = max(normwf);
[PPloc,NrPeaks,Ptail] = deal(nan(1,nwf_all));


for i = 1:nwf_all
    %Local PP is the peakiness of a part of the waveform; here, the part
    %surrounding the maximum peak. Local PP is used to distinguish between true
    %specular waveforms and specular waveforms with a tail
    if IDXpeak(i) > width && IDXpeak(i) < NrBins-width
        PPloc(i) = sum(wf_all(IDXpeak(i)-width:IDXpeak(i)+width,i))/sum(wf_all(:,i));
    else
        PPloc(i) = 10;
    end
    
    %Assess whether there are more significant peaks beside the main peak. The
    %significance of peaks is assessed based on their relative height and
    %prominence, and distance to the main peak. Peaks that are further away
    %have to be less high.
    NrPeaks(i) = numel(findpeaks(max(0,filter(-smooth_diff(8),1,wf_all(:,i))),1:NrBins,'MinPeakProminence',0.01,'MinPeakDistance',10));
    
    %Capture additional peaks preceding the main peak.
    if NrPeaks(i) == 1 && IDXpeak(i)<126
        NrPeaks(i) = numel(findpeaks(wf_all(1:IDXpeak(i)+2,i),1:IDXpeak(i)+2,'MinPeakHeight',0.5,'MinPeakProminence',0.05)); 
    end

    %Sum of total power in the tail (bins after the main peak+3). The value 
    %has to be neglectable for lead waveforms, but significant for ocean
    %ones. 
    Ptail(i)   = sum(wf_all(IDXpeak(i)+3:end,i));
    
end



%


%PPL & PPR 
for i = 1:nwf_all
    [~,locMAX] = max(normwf(:,i));
    PPL_all(i) = max(normwf(:,i))*3/(normwf(locMAX-3,i)+normwf(locMAX-2,i)+normwf(locMAX-1,i)) ; 
    PPR_all(i) = max(normwf(:,i))*3/(normwf(locMAX+3,i)+normwf(locMAX+2,i)+normwf(locMAX-+1,i)) ; 
end

%Waveform width 
for i = 1:nwf_all
    ww_all(i) = length(find (normwf(:,i)> 0.01*max(normwf(:,i)))) ; 
end
%Leading edge width
for i = 1:nwf_all
    MAX_12_bin_list = find(normwf(:,i)> 0.125*MAX_all(i)) ; 
    LeS_all(i) = MAX_bin(i) - MAX_12_bin_list(1) ; 
end
%TeS
for i = 1:nwf_all
    MIN_12_bin_list = find(normwf(:,i)< 0.125*MAX_all(i)) ; 
    TeS_all(i) = MAX_bin(i) - min(MIN_12_bin_list(MIN_12_bin_list > MAX_bin(i))); 
end

leadwf = wf_all(:,conf_all ==3) ; 
icewf = wf_all(:,conf_all ==1 ) ; 
maxleadwf = max(leadwf) ; 

figure 
subplot(2,3,1)
plot(leadwf(:,(maxleadwf < 1000)))
subplot(2,3,2)
plot(leadwf(:,(maxleadwf > 1000 & maxleadwf <5000)))
subplot(2,3,3)
plot(leadwf(:,(maxleadwf > 5000 & maxleadwf < 10000)))
subplot(2,3,4)
plot(leadwf(:,(maxleadwf > 10000 & maxleadwf <50000)))
subplot(2,3,5)
plot(leadwf(:,(maxleadwf > 50000 & maxleadwf <100000)))
subplot(2,3,6)
plot(leadwf(:,(maxleadwf > 100000 )))

% Normalize all the parameters ! 
kurt_all = kurt_all/max(kurt_all) ; 
skew_all = skew_all/max(skew_all) ; 
ww_all = ww_all/max(ww_all) ; 
LeS_all = LeS_all/max(LeS_all);
TeS_all = TeS_all/max(TeS_all);
sigma0_all = sigma0_all/max(sigma0_all) ; 


nwf_test = ceil(nwf_all * 0.20) ; 
idx_test   = randperm(nwf_all, nwf_test);
wf_test = wf_all(idx_test);
MAX_test = MAX_all(idx_test) ; 
kurt_test = kurt_all(idx_test);
skew_test = skew_all(idx_test) ; 
PP_test = PP_all(idx_test);
ww_test = ww_all(idx_test);
LeS_test = LeS_all(idx_test);
TeS_test = TeS_all(idx_test);
sigma0_test = sigma0_all(idx_test);
conf_test = conf_all(idx_test) ; 
PPL_test = PPL_all(idx_test) ; 
PPR_test = PPR_all(idx_test) ; 

testing_model = [MAX_test; kurt_test; PP_test ; ww_test ; LeS_test ; TeS_test;sigma0_test;  PPL_test; PPR_test].' ; 

idx_train = setdiff(1:nwf_all , idx_test) ; 
wf_train = wf_all(idx_train);
MAX_train = MAX_all(idx_train) ; 
kurt_train = kurt_all(idx_train);
skew_train = skew_all(idx_train) ; 
PP_train = PP_all(idx_train);
ww_train = ww_all(idx_train);
LeS_train = LeS_all(idx_train);
TeS_train = TeS_all(idx_train);
sigma0_train = sigma0_all(idx_train);
conf_train = conf_all(idx_train) ; 
PPL_train = PPL_all(idx_train) ; 
PPR_train = PPR_all(idx_train) ; 

training_model = [MAX_train; kurt_train; PP_train ; ww_train ; LeS_train ; TeS_train ; sigma0_train; PPL_train; PPR_train; conf_train].' ; 

train_all = [MAX_all; kurt_all; PP_all ; ww_all ; LeS_all ; TeS_all ; sigma0_all; PPL_all; PPR_all; PPloc; Ptail; NrPeaks; conf_all].' ; 
%%
% 
% 
% predictions = trainedModel.predictFcn(testing_model).' ; 
% iscorrect = predictions == conf_test ; 
% accuracy = sum(iscorrect)*100 / nwf_test ; %this is the overall accuracy for the test data set 
% % now find the TLR (which is actually the most important for this study !)
% % find true ice/lead & false ice/lead 
% 
% TrueIce = sum(predictions == 1 & conf_test == 1) ; 
% TrueLead = sum(predictions == 3 & conf_test == 3) ; 
% FalseIce = sum(predictions == 1 & conf_test == 3) ; 
% FalseLead = sum(predictions == 3 & conf_test == 1) ; 
% 
% TLR = TrueLead / (TrueLead + FalseIce) ; 
% FLR = FalseLead / (FalseLead + TrueIce) ; 
% % 
