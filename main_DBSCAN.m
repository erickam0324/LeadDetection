%% Script that uses DCSBAN as the classifier 
%% Need to merge files later so that the classifier and distance computation can be selected !!!!!! DAYO 
global D 
%% READ S3 L1b data
clear all; 
LoadCommonSettings_ericka;

rng(19970324)

SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          

FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';
[~,CS] = S3_L1b_read_ericka(FName)   ;
%% Code to raed OLCI GEO TIFF IMAGE 

% Locate OLCI Geo Tiff image : 
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
OLCIfile = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_RGB.tif") ;

GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
% OLCI = GEOTIFF_READ(OLCIfile);
% OLCILATmin = min(OLCI.y);
% OLCILATmax = max(OLCI.y);
% OLCILONmin = min(OLCI.x) ; 
% OLCILONmax = max(OLCI.x) ; 

Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
lat = ncread(GEOfile, "latitude");
lon = ncread(GEOfile, "longitude");

for i = 1:length(lon(:,1))
    for j = 1:length(lon(1,:))
        if lon(i,j) < 0 
            lon(i,j) = lon(i,j) + 360;
        end 
    end
end

% Find index of SAR data which corresponds to the image coordinates 
LONindex = find(CS.GEO.LON > min(lon(:)) & CS.GEO.LON < max(lon(:)) );
LATindex =  find(CS.GEO.LAT > min(lat(:)) & CS.GEO.LAT < max(lat(:)) ); 
validationindex = intersect(LONindex, LATindex) ; 


%% Classify waveforms 
DDAcf   = DDA_ConfigFile(SAT,'SAR');
Wf_Norm_Aw = 1;      
NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor.* CS.SAR.data;
unNORMdata =  CS.SAR.data;

% Apply SAR L1b_to_L2 code to find sigma0 
%[DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi) ; 
sigma0 =  CS.sigma0.' ; 

%waveforms that we are interested 
waveform = CS.SAR.data(:,validationindex) ; 
epsilon = 0.01; 
MinPts = 5 ; 
[class, isnoise, D] = DBSCAN(waveform,epsilon,MinPts); 
nclass = length(unique(class)); % number of class that was found by DBSCAN

figure 
for i = 0:nclass-1
     subplot(ceil(nclass/5),5,i+1)
     ind = find(class==i); 
     plot(waveform(:,ind))
     hold on
end
