%% Validation 

% algo flow : 
% select study area from validationoutput file 
% extract corresponding waveform from that study area 
% conduct classificatin 
% Compare results with validation (only with confident data points)


clear all; 
close all;
rng(19970324)

% retrieve all the necessary info 
D = '/Users/ericka/Desktop/Thesis/ValidationImageOutput';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    LeadConfidence{ii} = load(fullfile(D,N{ii},'lead_confidence.mat')); % improve by specifying the file extension. 
    wfSA{ii} = load(fullfile(D,N{ii},'wfSA.mat'));
    sigma0{ii} = load(fullfile(D,N{ii},'sigma0.mat'));
    SAR_idx{ii} = load(fullfile(D,N{ii},'SAR_idx.mat'));
    latSAR{ii} = load(fullfile(D,N{ii},'latSAR.mat'));
    lonSAR{ii} = load(fullfile(D,N{ii},'lonSAR.mat')) ; 
end


SA = 16 ; 
valFilename = N{SA}; 
OLCIsearchTag = valFilename(1:15); 

    %number of the study area that will be analyzed 
    predictions = fun_Kmed_KNN_classification(wfSA{SA}.wftosave, 30) ; 
    LeadConfidence{SA}.lead_confidence(LeadConfidence{SA}.lead_confidence==4) = 3 ; 

    %% ignore non-confident class 
    conf = LeadConfidence{SA}.lead_confidence ; 
    nonconf = conf(LeadConfidence{SA}.lead_confidence == 2) ; 
    predictions_conf = predictions ; 
    predictions_conf(LeadConfidence{SA}.lead_confidence == 2) = nan ; 

    % now find the TLR (which is actually the most important for this study !)
    % find true ice/lead & false ice/lead 

    TrueIce = sum(predictions_conf == 1 & conf == 1) ; 
    TrueLead = sum(predictions_conf == 3 & conf == 3) ; 
    FalseIce = sum(predictions_conf == 1 & conf == 3) ; 
    FalseLead = sum(predictions_conf == 3 & conf == 1) ; 
    accuracy = (TrueIce + TrueLead) / (TrueIce + TrueLead + FalseIce + FalseLead) ; 

    TLR = TrueLead / (TrueLead + FalseIce) ; 
    FLR = FalseLead / (FalseLead + TrueIce) ; 
    % 
    
    


    %% visualize 


        OLCIpath = '/Users/ericka/Desktop/Thesis/OLCIdata';
        S1 = dir(fullfile(OLCIpath,'*'));
        NamesOLCIdata = setdiff({S1([S1.isdir]).name},{'.','..'}); % list of subfolders of D.
        SEN3file = NamesOLCIdata{contains(NamesOLCIdata, OLCIsearchTag)==1 } ; 
        %SEN3file = 'S3A_OL_1_EFR____20170312T010959_20170312T011201_20180410T134742_0122_015_188_1620_LR2_R_NT_002.SEN3'; 
        date = SEN3file(17:24) ; 
        SEN3file = fullfile(OLCIpath, SEN3file) ;

        Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
        Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
        Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");

        GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
        lat = ncread(GEOfile, "latitude");
        lon = ncread(GEOfile, "longitude");

        A_col(:,:,1) = Oa08;  % red
        A_col(:,:,2) = Oa05;  % green
        A_col(:,:,3) = Oa03;  % blue


    figure, worldmap([min(latSAR{SA}.latSARtosave) max(latSAR{SA}.latSARtosave)], [min(lonSAR{SA}.lonSARtosave) max(lonSAR{SA}.lonSARtosave)])
    geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
    geoimg.AlphaDataMapping = 'none';
    geoimg.FaceAlpha = 'texturemap';
    alpha(geoimg,double(~isnan(Oa03)))
    hold on 
    % Truelead = R, Trueice = B 
    % FalseLead = m, FalseIce = c
    scatterm(latSAR{SA}.latSARtosave(predictions_conf == 1 & conf == 1),lonSAR{SA}.lonSARtosave(predictions_conf == 1 & conf == 1), 'filled', 'b') ;
    scatterm(latSAR{SA}.latSARtosave(predictions_conf == 1 & conf == 3),lonSAR{SA}.lonSARtosave(predictions_conf == 1 & conf == 3), 'filled', 'c') ;
    scatterm(latSAR{SA}.latSARtosave(predictions_conf == 3 & conf == 1),lonSAR{SA}.lonSARtosave(predictions_conf == 3 & conf == 1), 'filled', 'm') ;
    scatterm(latSAR{SA}.latSARtosave(predictions_conf == 3 & conf == 3),lonSAR{SA}.lonSARtosave(predictions_conf == 3 & conf == 3), 'filled', 'r') ;

    
    %% to change index of the validation labels : uncomment the following and run this 
    
    
% save(fullfile(D,N{SA},'lead_confidencenew.mat'),'lead_confidence')
    