%% Trying to plot a portion of OLCI image (only the relevant areas)
clear all; 
LoadCommonSettings_ericka;
rng(19970324)
%% Load OLCI related files
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance");
Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
Oa15 = ncread(fullfile(SEN3file, "Oa15_radiance.nc"), "Oa15_radiance");
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

FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';
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
lat_max_des = 77.55;
lat_min_des = 77.25;
lon_max_des = 180.0;
lon_min_des = 179.5;
crop_idx_lat = find(lat < lat_min_des | lat > lat_max_des ) ; 
crop_idx_lon = find(lon < lon_min_des | lon > lon_max_des ) ;
crop_idx = intersect(crop_idx_lat, crop_idx_lon)  ;  % these are the indices that should be cropped !

Oa08(crop_idx) = NaN; 
Oa05(crop_idx) = NaN;
Oa03(crop_idx) = NaN;
Oa01(crop_idx) = NaN;

NDVI = (Oa01-Oa21)./(Oa21+Oa01); 
NDVI = NDVI./(max(NDVI(:))); 

% Plot
% figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
%figure, worldmap([min(LAT_SAR) max(LAT_SAR)], [min(LON_SAR) max(LON_SAR)])

A(:,:,1) = Oa03;  % red
A(:,:,2) = Oa05;  % green
A(:,:,3) = Oa08;  % blue

tic 
im = A; 
imflat = double(reshape(normalize(im), size(im,1) * size(im,2),3)) ; 
K = 2; 
[kIDs, kC] = kmeans(imflat, K, 'Display', 'iter', 'MaxIter', 150, 'Start', 'sample') ;
imout = reshape(kIDs, size(im,1), size(im,2)) ;
toc 

imout(imout==1) = 0;
imout(imout==2) = 60;
imout(imout==3) = 150;
imout(imout==4) = 200;
imout(imout==5) = 250;


%figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
% figure, worldmap([76.5 77.8], [176.5 180.6])
figure, worldmap([lat_min_des lat_max_des], [lon_min_des lon_max_des])
geoimg=geoshow(lat, lon, imout/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(CS.GEO.LAT(OLCI_intersec_idx(iceindex)),CS.GEO.LON(OLCI_intersec_idx(iceindex)), 'filled', 'b'); 
scatterm(CS.GEO.LAT(OLCI_intersec_idx(leadindex)),CS.GEO.LON(OLCI_intersec_idx(leadindex)),  'filled', 'r'); 


