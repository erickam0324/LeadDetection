%% loading icebrige tiff file 
%I = GEOTIFF_READ('icebridge/122690515/DMS_1742629_08917_20170419_15544407.tif');

[a,R] = geotiffread('icebridge/122690515/DMS_1742629_08917_20170419_15544407.tif');
figure
mapshow(a,R);
type('icebridge/122690515/DMS_1742629_08917_20170419_15544407.tif.xml')

lon_ib1 = [-122.357783  -122.351771 -122.322520 -122.328545]; 
lat_ib1 = [77.541628  77.547410  77.545992  77.540211] ; 

lon_ib2 = [-122.541215  -122.535243 -122.505843 -122.511828]; 
lat_ib2 = [77.567318  77.573137  77.571736  77.565918] ; 

lon_ib3 = [-122.037501  -122.031424 -122.002544 -122.008634]; 
lat_ib3 = [77.497393  77.503108  77.501667  77.495952] ; 

lon_ib4 = [-122.106499  -122.100415 -122.071401 -122.077498]; 
lat_ib4 = [77.506535  77.512284  77.510845  77.505096] ; 


%% OLCI 
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = fullfile(OLCIpath, "S3A_OL_1_EFR____20170419T195450_20170419T195750_20180416T235955_0179_016_356_1620_LR2_R_NT_002.SEN3") ;
GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance") ;
Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance") ;
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance") ;
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance") ;
lat = ncread(GEOfile, "latitude") ;
lon = ncread(GEOfile, "longitude") ;

% Make negative longitude values positive and express angles between 0 and
% 360
for i = 1:length(lon(:,1))
    for j = 1:length(lon(1,:))
        if lon(i,j) < 0 
            lon(i,j) = lon(i,j) + 360;
        end 
    end
end

A(:,:,1) = Oa08;  % red
A(:,:,2) = Oa05;  % green
A(:,:,3) = Oa03;  % blue

figure, worldmap([77 78],[235 239]); 
geoimg=geoshow(lat, lon, A/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(lat_ib1,lon_ib1,'filled', 'r')
scatterm(lat_ib2,lon_ib2,'filled', 'k')
scatterm(lat_ib3,lon_ib3,'filled', 'b')
scatterm(lat_ib4,lon_ib4,'filled', 'g')



