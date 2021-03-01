%% Threshold overlaying greyscale 

%% READ S3 L1b data
clear all; 
LoadCommonSettings_ericka;

SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          

FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';
[~,CS] = S3_L1b_read_ericka(FName)   ;



%% 
%% Code to raed OLCI GEO TIFF IMAGE 
% Locate original olci image. 
% The coordinates are going to be the same as the nc files. 

OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3") ;
%OLCIfile = fullfile(OLCIpath, "S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_Oa09_radiance.tif");

GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
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
NORMfactor = NORMfactor(validationindex) ; 
unNORMdata =  CS.SAR.data(:,validationindex);

% Apply SAR L1b_to_L2 code to find sigma0 
[DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi) ; 
sigma0 =  CS.sigma0(validationindex).' ; 

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
%% Image segmentation 
imout = OLCIsegmentation('/Users/ericka/Desktop/Thesis/OLCIData/olci.png', "graytif") ; 

imout = OLCIsegmentation('/Users/ericka/Desktop/Thesis/OLCIData/S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_Oa10_radiance.png', "graytif") ; 
% first input is the image retrieved from the snap software 
% second input is the desired name for the new grayscale tiff image 

%% Plot new grayscale in figure 
tiffile = '/Users/ericka/Desktop/Thesis/OLCIData/S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002_Oa08_radiance.tif'; 
[~,R,~] = readgeoraster(tiffile) ; 
R.LatitudeLimits = [min(lat(:)), max(lat(:))] ; 
R.LongitudeLimits = [min(lon(:)), max(lon(:))] ; 
geotiffwrite('/Users/ericka/Desktop/Thesis/output/greyscaletif/graytif.tif', imout-1, colormap(bone(3)), R)


figure, worldmap([min(lat(:)) max(lat(:))],[min(lon(:)) max(lon(:))])
geoshow(lat,lon,(imout-1).', colormap(bone(3)))
hold on 
scatterm(CS.GEO.LAT(validationindex(iceindex)),CS.GEO.LON(validationindex(iceindex)), 'filled', 'b'); 
scatterm(CS.GEO.LAT(validationindex(leadindex)),CS.GEO.LON(validationindex(leadindex)),  'filled', 'r'); 

 

