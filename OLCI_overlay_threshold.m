%% READ S3 L1b data
clear all; 
LoadCommonSettings_ericka;


SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          

SARpath = '/Users/ericka/Desktop/Thesis/SARdata';
FName = 'S3A_SR_1_SRA____20170312T010714_20170312T015743_20170406T155339_3029_015_188______MAR_O_NT_002.SEN3';
FName = fullfile(SARpath, FName) ;
[~,CS, NORMfactor] = S3_L1b_read_ericka(FName)   ;



%% Code to raed OLCI GEO TIFF IMAGE 

% Locate OLCI Geo Tiff image : 
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
%SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
SEN3file = fullfile(OLCIpath, 'S3A_OL_1_EFR____20170312T010959_20170312T011201_20180410T134742_0122_015_188_1620_LR2_R_NT_002.SEN3') ;

GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 


Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance");
Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
Oa21 = ncread(fullfile(SEN3file, "Oa21_radiance.nc"), "Oa21_radiance");
Oa15 = ncread(fullfile(SEN3file, "Oa15_radiance.nc"), "Oa15_radiance");
lat = ncread(GEOfile, "latitude");
lon = ncread(GEOfile, "longitude");

% LOOKING AT WHERE OLCI IMAGE IS LOCATED. 
% if olci is has coordinates that are crossing values at 0deg, the
% coordinates should be expressed in -180 to 180. Otherwise, 0 to 360 
if ~isempty(find(lon>0 & lon<1, 1)) && ~isempty(find(lon<0 & lon> -1, 1))
    % use -180 to 180 coordinate system (change SAR lon values)
    coordinate = '-180to180' ; 
        for i = 1:length(CS.GEO.LON)
            if CS.GEO.LON(i) >180
               CS.GEO.LON(i) =  CS.GEO.LON(i) - 360;
            end 
        end
else 
    % use 0 to 360 coordinate system (change OLCI lon values)
    coordinate = '0to360';
        for i = 1:length(lon(:,1))
            for j = 1:length(lon(1,:))
                if lon(i,j) < 0 
                    lon(i,j) = lon(i,j) + 360;
                end 
            end
        end
end 



% Find index of SAR data which corresponds to the image coordinates 
LONindex = find(CS.GEO.LON > min(lon(:)) & CS.GEO.LON < max(lon(:)) );
LATindex =  find(CS.GEO.LAT > min(lat(:)) & CS.GEO.LAT < max(lat(:)) ); 
OLCI_intersec_idx = intersect(LONindex, LATindex) ; 



%% Classify waveforms 
NORMfactor = NORMfactor(OLCI_intersec_idx) ; 
unNORMdata =  CS.SAR.data(:,OLCI_intersec_idx);
sigma0 =  CS.sigma0(OLCI_intersec_idx).' ; 

classification_method = 'Inger'; % choose from 'Inger', 'Rose', 'Peacock', ... 
class = Classify_Waveform(unNORMdata, NORMfactor, sigma0, classification_method) ;  

%% 
leadindex =  find(class == 2 ); 
iceindex = find(class == 1) ; 
ambindex = find(class == 0 ); 
% % 
%   geoscatter(CS.GEO.LAT(iceindex),CS.GEO.LON(iceindex), 'filled', 'b')
%   hold on
%   geoscatter(CS.GEO.LAT(leadindex),CS.GEO.LON(leadindex),  'filled', 'r')
%   %scatter( CS.GEO.LAT(ambindex),CS.GEO.LON(ambindex), 'filled', 'g')
%   hold off
 
%%
% create true or false 3-color image
clear A ; 
A(:,:,1) = Oa08;  % red
A(:,:,2) = Oa05;  % green
A(:,:,3) = Oa03;  % blue

NDVI = (Oa01-Oa21)./(Oa21+Oa01); 
NDVI = NDVI./(max(NDVI(:))); 
%%

figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
%figure, worldmap([77 83], [330 360])
geoimg=geoshow(lat, lon, A/256,'DisplayType','image');  
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(CS.GEO.LAT(OLCI_intersec_idx(iceindex)),CS.GEO.LON(OLCI_intersec_idx(iceindex)), 'filled', 'b'); 
scatterm(CS.GEO.LAT(OLCI_intersec_idx(leadindex)),CS.GEO.LON(OLCI_intersec_idx(leadindex)),  'filled', 'r'); 
title(date)
 
