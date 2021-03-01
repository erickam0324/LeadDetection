%% ValidationScheme 

%% Trying to plot a portion of OLCI image (only the relevant areas)
clear all; 
LoadCommonSettings_ericka;
rng(19970324)
%% Load OLCI related files
OLCIpath = "/Users/ericka/Desktop/Thesis/OLCIData" ; 
SEN3file = 'S3A_OL_1_EFR____20170401T005019_20170401T005319_20180413T161703_0179_016_088_1620_LR2_R_NT_002.SEN3' ; 
date = SEN3file(17:24) ; 
SEN3file = fullfile(OLCIpath, SEN3file) ;

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
% enter following lon MAX MIN showing the study area. These values may change
% slightly depending on the sizes of the small squares that will be used 
% for image segmentation analysis. 
lon_max_SA = 180.5;
lon_min_SA = 177.7;

% count how many squares we need to analyse 
lon_square = 0.6 ; % size (degree of longitude) of a single square
patch_size = 0.1 ; % this is the size (degree of longitude) that two squares overlap such that validation scheme can be carried out for all points (if 0. the edge points might not have enough data point for validation)
N_square = ceil((lon_max_SA - lon_min_SA) / (lon_square-patch_size)) ; 

% Make list of max/min lon/lat per square with overlaps 
for i = 1:N_square
   lon_max_list(i) = lon_max_SA - (lon_square-patch_size)*(i-1); 
   lon_min_list(i) = lon_max_list(i) - lon_square ; 
   lat_max_list(i) = CS.GEO.LAT(find(CS.GEO.LON < lon_max_list(i), 1 )) + 0.05; 
   lat_min_list(i) = CS.GEO.LAT(find(CS.GEO.LON > lon_min_list(i), 1, 'last' )) - 0.05;
end 

lat_max_SA = max(lat_max_list);
lat_min_SA = min(lat_min_list);
lon_max_SA = max(lon_max_list);
lon_min_SA = min(lon_min_list);

A_col(:,:,1) = Oa08;  % red
A_col(:,:,2) = Oa05;  % green
A_col(:,:,3) = Oa03;  % blue

%% image segmentation per squaere(Kmeans)

for i = 1:N_square
    % find indices that should be cropped per square for faster computation
    crop_idx_lat_SA = find(lat < lat_min_list(i) | lat > lat_max_list(i) ) ; 
    crop_idx_lon_SA = find(lon < lon_min_list(i) | lon > lon_max_list(i) ) ;
    crop_idx_SA = intersect(crop_idx_lat_SA, crop_idx_lon_SA)  ;  
    
    Oa08_SA = Oa08;
    Oa05_SA = Oa05;
    Oa03_SA = Oa03; 
    Oa15_SA = Oa15; 
    Oa21_SA = Oa21;
    Oa08_SA(crop_idx_SA) = NaN; 
    Oa05_SA(crop_idx_SA) = NaN;
    Oa03_SA(crop_idx_SA) = NaN;
    Oa15_SA(crop_idx_SA) = NaN;
    Oa21_SA(crop_idx_SA) = NaN;
    
    A(:,:,1) = Oa03_SA;  % red
    A(:,:,2) = Oa05_SA;  % green
    A(:,:,3) = Oa08_SA;  % blue
    
    % Kmeans segmentation 
    tic 
    im = A; 
    imflat = double(reshape(normalize(im), size(im,1) * size(im,2),3)) ; 
    K = 2; 
    [kIDs, kC] = kmeans(imflat, K, 'Display', 'iter', 'MaxIter', 150, 'Start', 'sample') ;
    imout = reshape(kIDs, size(im,1), size(im,2)) ;
    toc 
    
    % Make the colors consistent. Assuming that there are less 'lead pixel'
    % in an image compared to 'ice pixel', the image will show darker color
    % for minority pixels : 
    
    if length(find(imout==1)) >= length(find(imout==2))
        ind_1 = find(imout==1) ; 
        ind_2 = find(imout==2) ; 
        imout(ind_1) = 2; 
        imout(ind_2) = 1; 
        clear ind_1 ind_2
    end 
    
    imout(imout==1) = 0; % Leads, black
    imout(imout==2) = 50; % Ice, grey 
    
    % find N closest pixels to the SAR track from OLCI image (within the current square)
    Nclose_pixel = 3; 
    [idx_OLCIunderSAR{i}, SAR_idx_persquare{i}] = correspondingOLCI(CS, lat, lon, lon_max_list(i), lon_min_list(i), Nclose_pixel ) ; 
    
    % find whether these N closest pixels are leads or ice 
    for j = 1:length(idx_OLCIunderSAR{i})
        result_N{i}{j} = imout([idx_OLCIunderSAR{i}{j}]) ; 
        if length(find(result_N{i}{j}==0)) >= ceil(Nclose_pixel/2)
            result_ave{i}{j} = 2 ; % LEAD = 2
        elseif length(find(result_N{i}{j}==50)) >= ceil(Nclose_pixel/2)
            result_ave{i}{j} = 1 ; % ICE = 1
        else 
            result_ave{i}{j} = 0 ; % AMBIGUOUS (something is wrong ... )
        end 
    end 
    
    
    % plot the segmented image on the left 
%     subplot(1,2,1), worldmap([lat_min_list(i) lat_max_list(i)], [lon_min_list(i) lon_max_list(i)])
%     geoimg=geoshow(lat, lon, imout/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
%     geoimg.AlphaDataMapping = 'none';
%     geoimg.FaceAlpha = 'texturemap';
%     alpha(geoimg,double(~isnan(Oa03)))
%     hold on 
%     scatterm(CS.GEO.LAT(OLCI_intersec_idx(iceindex)),CS.GEO.LON(OLCI_intersec_idx(iceindex)), 'filled', 'b'); 
%     scatterm(CS.GEO.LAT(OLCI_intersec_idx(leadindex)),CS.GEO.LON(OLCI_intersec_idx(leadindex)),  'filled', 'r'); 
%     
%     % plot the original image on the right
%     subplot(1,2,2), worldmap([lat_min_list(i) lat_max_list(i)], [lon_min_list(i) lon_max_list(i)])
%     geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  % this 200 was a bit random, depends on the maximum values in the bands
%     geoimg.AlphaDataMapping = 'none';
%     geoimg.FaceAlpha = 'texturemap';
%     alpha(geoimg,double(~isnan(Oa03)))
%     hold on 
%     scatterm(CS.GEO.LAT(OLCI_intersec_idx(iceindex)),CS.GEO.LON(OLCI_intersec_idx(iceindex)), 'filled', 'b'); 
%     scatterm(CS.GEO.LAT(OLCI_intersec_idx(leadindex)),CS.GEO.LON(OLCI_intersec_idx(leadindex)),  'filled', 'r'); 
%     
%     % save the image 
%     ROOT_imageoutput = '/Users/ericka/Desktop/Thesis/ValidationImageOutput' ; 
%     if not(exist([ROOT_imageoutput,'/',date,'_LON',num2str(lon_max_SA),'_',num2str(lon_min_SA)]))
%         mkdir(['/Users/ericka/Desktop/Thesis/ValidationImageOutput/',date,'_LON',num2str(lon_max_SA),'_',num2str(lon_min_SA)])
%     end 
%     saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_max_SA),'_',num2str(lon_min_SA),'/image',num2str(i),'_','LON',num2str(lon_min_list(i)),'_',num2str(lon_max_list(i)),'.png'])
%     
%     close all 
%     clear Oa08_SA Oa05_SA Oa03_SA Oa15_SA Oa21_SA

end

% we now have SAR_idx_persquare that shows the indicies of SAR of Studyarea
% and results_ave that shows the class of those points
% there's some over-lap but that will be resolved when making into one
% array and the results from the 'next' square will overwrite the result of
% the 'previous' square. 

% make the indicies in one array 
for i = 1:N_square 
    add1 = SAR_idx_persquare{i} ; 
    add2 = cell2mat(result_ave{i}).' ; 
    if i ==1 
        SAR_idx_all = [add1] ; 
        result_all = [add2] ; 
    else 
        SAR_idx_all = [SAR_idx_all; add1] ; 
        result_all = [result_all; add2] ; 

    end 
end 


    % plot
    figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
    geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
    geoimg.AlphaDataMapping = 'none';
    geoimg.FaceAlpha = 'texturemap';
    alpha(geoimg,double(~isnan(Oa03)))
    hold on 
    scatterm(CS.GEO.LAT(SAR_idx_all(result_all==1)),CS.GEO.LON(SAR_idx_all(result_all==1)), 'filled', 'b'); 
    scatterm(CS.GEO.LAT(SAR_idx_all(result_all==2)),CS.GEO.LON(SAR_idx_all(result_all==2)),  'filled', 'r'); 
    
%% Removing obvious inaccurate points 
    % now there might be points that are classified into a class that are
    % clearly wrong (refer to the segmentation images)
    % a prompt asking whether there should be a part that should be 'unclassified'

prompt_YN = 'Are there areas (still) that should be discarded? Y/N : ';
YN = input(prompt_YN, 's') ; 
i = 1; 
while YN =='Y'
   prompt_discardlon = 'Enter range of longitude that should be discarded [min_lon, max_lon] (eg: [172.3, 172.6]): ';
   discard_lon(i,:) = input(prompt_discardlon) ; 
   i = i +1; 
   YN = input(prompt_YN, 's') ;
end 
 

% find indicies of SAR altimetry that should be discarded 
for i = 1: length(discard_lon(:,1))
    index_to_discard(i,:) = find(CS.GEO.LON > discard_lon(i,1) & CS.GEO.LON < discard_lon(i,2) ) ;    
end 

% discard (=nan) the results that corresponds to those indicies 
result_all(ismember(SAR_idx_all ,index_to_discard)==1) = nan ; 

%% 
% dealing with the overlapps
% making sure that the repeated indicies are the ones that the results will
% be considered. (indicies that appear for the second time)

[SAR_idx_all_nooverlap, ~, ic] = unique(SAR_idx_all, 'last') ; 
result_all_nooverlap(ic) = result_all ;



