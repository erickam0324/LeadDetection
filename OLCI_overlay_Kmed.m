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


SARpath = "/Users/ericka/Desktop/Thesis/SARdata" ; 
FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';
%FName = 'S3B_SR_1_SRA____20200331T134818_20200331T143848_20200425T213147_3029_037_167______LN3_O_NT_003.SEN3' ; 
SEN3file = fullfile(SARpath, FName) ;
[~,CS] = S3_L1b_read_ericka(FName)   ;

%% Code to raed OLCI GEO TIFF IMAGE 

% Locate OLCI Geo Tiff image : 
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
%SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20200317T103030_20200317T103313_20200318T150854_0163_056_108_1620_LN1_O_NT_002.SEN3") ;
OLCIfile = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_RGB.tif") ;
%OLCIfile = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_Oa08_radiance.tif");
%OLCIfile = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_Oa09_radiance.tif");

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
NORMdata   = NORMfactor .* CS.SAR.data;
unNORMdata =  CS.SAR.data;

% Apply SAR L1b_to_L2 code to find sigma0 
[DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi) ; 
sigma0 =  CS.sigma0.' ; 

%waveforms that we are interested 
waveform = CS.SAR.data(:,validationindex) ; 
class = fun_Kmed_KNN_classification(waveform) ;  


leadindex =  find(class == 2 );
iceindex = find(class == 1 ) ;
ambindex = find(class == 0 ) ;
oceanindex = find(class == 3) ;
outlierindex = find(class == 4) ;

%%
% create true or false 3-color image
clear A ; 
A(:,:,1) = Oa08;  % red
A(:,:,2) = Oa05;  % green
A(:,:,3) = Oa03;  % blue

figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
%figure, worldmap([76.4 78.3], [176 184])
geoimg=geoshow(lat, lon, A/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(CS.GEO.LAT(validationindex(iceindex)),CS.GEO.LON(validationindex(iceindex)), 'filled', 'b'); 
scatterm(CS.GEO.LAT(validationindex(leadindex)),CS.GEO.LON(validationindex(leadindex)),  'filled', 'r'); 
scatterm(CS.GEO.LAT(validationindex(oceanindex)),CS.GEO.LON(validationindex(oceanindex)),  'filled', 'c');
% scatterm(CS.GEO.LAT(validationindex(ambindex)),CS.GEO.LON(validationindex(ambindex)),  'filled', 'k'); 
% scatterm(CS.GEO.LAT(validationindex(ambindex)),CS.GEO.LON(validationindex(ambindex)),  'filled', 'g'); 

%% 
figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
geoimg=geoshow(lat, lon, Oa08,'DisplayType','image');

