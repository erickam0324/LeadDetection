%% Validation scheme 2 
%% looking at sudden changes of radiance along track of altimetry measurements 
%% Trying to plot a portion of OLCI image (only the relevant areas)
clear all; 
LoadCommonSettings_ericka;
rng(19970324)
%% Load OLCI related files
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
%Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance");
Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
Oa15 = ncread(fullfile(SEN3file, "Oa15_radiance.nc"), "Oa15_radiance");
Oa20 = ncread(fullfile(SEN3file, "Oa20_radiance.nc"), "Oa20_radiance");
Oa21 = ncread(fullfile(SEN3file, "Oa21_radiance.nc"), "Oa21_radiance");
lat = ncread(GEOfile, "latitude");
lon = ncread(GEOfile, "longitude");

% Make negative longitude values positive and express angles between 0 and
% 360
for i = 1:length(lon(:,1))
    for j = 1:length(lon(1,:))
        if lon(i,j) < 0 
            lon(i,j) = lon(i,j) + 360;
        end 
    end
end

%% Load SAR altimetery files 
SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          


SARpath = "/Users/ericka/Desktop/Thesis/SARdata" ; 
FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';
%FName = 'S3B_SR_1_SRA____20200331T134818_20200331T143848_20200425T213147_3029_037_167______LN3_O_NT_003.SEN3' ; 
SEN3file = fullfile(SARpath, FName) ;
[~,CS] = S3_L1b_read_ericka(FName)   ;

% Find index of SAR data which corresponds to the image coordinates 
LONindex = find(CS.GEO.LON > min(lon(:)) & CS.GEO.LON < max(lon(:))) ;
LATindex =  find(CS.GEO.LAT > min(lat(:)) & CS.GEO.LAT < max(lat(:)) ); 
OLCI_intersec_idx = intersect(LONindex, LATindex) ; 
LON_SAR = CS.GEO.LON(OLCI_intersec_idx);
LAT_SAR = CS.GEO.LAT(OLCI_intersec_idx);

%% Classify waveforms 
DDAcf   = DDA_ConfigFile(SAT,'SAR');
Wf_Norm_Aw = 1;      
NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
NORMfactor = NORMfactor(OLCI_intersec_idx) ; 
unNORMdata =  CS.SAR.data(:,OLCI_intersec_idx);

% Apply SAR L1b_to_L2 code to find sigma0 
[DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi) ; 
sigma0 =  CS.sigma0(OLCI_intersec_idx).' ; 

classification_method = 'Inger'; % choose from 'Inger', 'Rose', 'Peacock', 'Kmed'
class = Classify_Waveform(unNORMdata, NORMfactor, sigma0, classification_method) ;  

leadindex =  find(class == 2 ); 
iceindex = find(class == 1) ; 
ambindex = find(class == 0 ); 

%% RUN IF WANT TO CROP Cropping out irrelevant indicies in the image matrix 
lon_max_SA = 180.5;
lon_min_SA = 177.4;

[lat_min_SA, lat_max_SA] = correspondingLAT(CS, lon_min_SA, lon_max_SA, 0.05) ; % last input = margin, how much margin you want the output map to have  
crop_idx_lat = find(lat < lat_min_SA | lat > lat_max_SA ) ; 
crop_idx_lon = find(lon < lon_min_SA | lon > lon_max_SA ) ;
crop_idx = intersect(crop_idx_lat, crop_idx_lon)  ;  % these are the indices that should be cropped !

Oa21(crop_idx) = NaN; 
Oa20(crop_idx) = NaN;
Oa08(crop_idx) = NaN; 
Oa05(crop_idx) = NaN;
Oa03(crop_idx) = NaN;
Oa01(crop_idx) = NaN;


A(:,:,1) = Oa03;  % red
A(:,:,2) = Oa05;  % green
A(:,:,3) = Oa08;  % blue

%%
rad = Oa08 ;  % Oa08 or the NDSIII
Nclose_pixel = 1; 

[idx_OLCIunderSAR, ~] = correspondingOLCI(CS, lat, lon, lon_max_SA, lon_min_SA, Nclose_pixel ) ; 

rad_series = rad([idx_OLCIunderSAR{:}]) ; 

figure
plot(1:length(rad_series), rad_series)
[max_pks, max_pks_idx] = findpeaks(rad_series) ; 
[min_pks, min_pks_idx] = findpeaks(-rad_series) ; 

%min peaks are the possible leads (if large enough compared to previous max) 
%for every min peak, the difference to the previous max peak will be
% checked. 

%%
for i = 1:length(min_pks) 
    if isempty(find(max_pks_idx < min_pks_idx(i), 1)) 
        previous_max_idx(i) = NaN ; 
        previous_max(i) = NaN ; 
    else
        previous_max_idx(i) = max_pks_idx(find(max_pks_idx < min_pks_idx(i), 1, 'last')) ; 
        previous_max(i) = max_pks(find(max_pks_idx < min_pks_idx(i), 1, 'last')) ; 
    end 
    
    if isempty(find(max_pks_idx > min_pks_idx(i), 1))  
        next_max(i) = NaN ; 
        next_max_idx(i) = NaN ; 
    else
        next_max_idx(i) = max_pks_idx(find(max_pks_idx > min_pks_idx(i), 1)); 
        next_max(i) = max_pks(find(max_pks_idx > min_pks_idx(i), 1)) ; 
    end 
    
end 
%%
% threshold is given such that if the dip is larger than that ratio
% compared to the previous max peak, then that index is considered as a
% lead
threshold = 0.015;  
lead_min_pks = find(abs(min_pks./previous_max) < (1- threshold)); 
lead_idx = [];
for i = 1:length(lead_min_pks)
    j = lead_min_pks(i) ; 
    append = previous_max_idx(j)+1:next_max_idx(j)-1 ; % maybe add if statement if the lead is very dark it can be wider.. 
    lead_idx = [lead_idx append] ; 
end 

OLCI_lead_lat = lat([idx_OLCIunderSAR{lead_idx}]) ; 
OLCI_lead_lon = lon([idx_OLCIunderSAR{lead_idx}]) ; 

figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
geoimg=geoshow(lat, lon, A/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(OLCI_lead_lat, OLCI_lead_lon,  'filled', 'r'); 
%scatterm(CS.GEO.LAT(LONindex_SA(i)),CS.GEO.LON(LONindex_SA(i)),  'filled', 'r'); 
% scatterm(CS.GEO.LAT(OLCI_intersec_idx(iceindex)),CS.GEO.LON(OLCI_intersec_idx(iceindex)), 'filled', 'b'); 
% scatterm(CS.GEO.LAT(OLCI_intersec_idx(leadindex)),CS.GEO.LON(OLCI_intersec_idx(leadindex)),  'filled', 'r'); 
