%% validation all together !! 
% combining validation 1 & 2 into one script. 

clear all; 
LoadCommonSettings_ericka;
rng(19970324)

ii = 1 ; 

    OLCIpath = '/Users/ericka/Desktop/Thesis/OLCIdata';
    S1 = dir(fullfile(OLCIpath,'*'));
    NamesOLCIdata = setdiff({S1([S1.isdir]).name},{'.','..'}); % list of subfolders of D.

    SARpath = '/Users/ericka/Desktop/Thesis/SARdata';
    S2 = dir(fullfile(SARpath,'*'));
    NamesSARdata = setdiff({S2([S2.isdir]).name},{'.','..'}); % list of subfolders of D.

    FName = NamesSARdata{ii} ; 
    FName = fullfile(SARpath, FName) ;
    [~,CS, NORMfactor] = S3_L1b_read_ericka(FName)   ;

    % Code to raed OLCI GEO TIFF IMAGE 
    % Locate OLCI Geo Tiff image : 
    SEN3file = NamesOLCIdata{ii} ; 
    date = SEN3file(17:24) ; 
    SEN3file = fullfile(OLCIpath, SEN3file) ;
 %%
ROOT_imageoutput = '/Users/ericka/Desktop/Thesis/ValidationImageOutput' ; 

GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance");
Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
Oa15 = ncread(fullfile(SEN3file, "Oa15_radiance.nc"), "Oa15_radiance");
Oa21 = ncread(fullfile(SEN3file, "Oa21_radiance.nc"), "Oa21_radiance");
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

A_col(:,:,1) = Oa08;  % red
A_col(:,:,2) = Oa05;  % green
A_col(:,:,3) = Oa03;  % blue

% Load SAR altimetery files 
SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          


% Find index of SAR data which corresponds to the image coordinates 
LONindex = find(CS.GEO.LON > min(lon(:)) & CS.GEO.LON < max(lon(:))) ;
LATindex =  find(CS.GEO.LAT > min(lat(:)) & CS.GEO.LAT < max(lat(:)) ); 
OLCI_intersec_idx = intersect(LONindex, LATindex) ; 
LON_SAR = CS.GEO.LON(OLCI_intersec_idx);
LAT_SAR = CS.GEO.LAT(OLCI_intersec_idx);


%% Classify waveforms 
sigma0 =  CS.sigma0(OLCI_intersec_idx).' ; 
NORMfactor = NORMfactor(OLCI_intersec_idx) ; 

classification_method = 'Inger'; % choose from 'Inger', 'Rose', 'Peacock', 'Kmed'
class = Classify_Waveform(CS.SAR.data(:,OLCI_intersec_idx), NORMfactor, sigma0, classification_method) ;  

leadindex =  find(class == 2 ); 
iceindex = find(class == 1) ; 
ambindex = find(class == 0 ); 


%% ---------------------------- Study area selection ---------------------------------
% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------
% 
% enter following lon MAX MIN showing the study area. These values may change
% slightly depending on the sizes of the small squares that will be used 
% for image segmentation analysis. 
lon_min_SA = 160; 
lon_max_SA = 165;
lat_min_SA = min(CS.GEO.LAT(find(CS.GEO.LON > lon_min_SA & CS.GEO.LON < lon_max_SA )))-0.05 ; 
lat_max_SA = max(CS.GEO.LAT(find(CS.GEO.LON > lon_min_SA & CS.GEO.LON < lon_max_SA )))+0.05 ; 




%% ---------------------------- ValidationScheme 1 -----------------------------------
% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------

% count how many squares we need to analyse 
lon_square = 0.6 ; % size (degree of longitude) of a single square
patch_size = 0.1 ; % this is the size (degree of longitude) that two squares overlap such that validation scheme can be carried out for all points (if 0. the edge points might not have enough data point for validation)
N_square = ceil((lon_max_SA - lon_min_SA) / (lon_square-patch_size)) ; 

% Make list of max/min lon/lat per square with overlaps 
for i = 1:N_square
   lon_max_list(i) = lon_max_SA - (lon_square-patch_size)*(i-1); 
   lon_min_list(i) = lon_max_list(i) - lon_square ; 
   lat_max_list(i) = max(CS.GEO.LAT(find(CS.GEO.LON > lon_min_list(i) & CS.GEO.LON < lon_max_list(i) ))) +0.05 ; 
   lat_min_list(i) = min(CS.GEO.LAT(find(CS.GEO.LON > lon_min_list(i) & CS.GEO.LON < lon_max_list(i) ))) -0.05 ; 
end 

% New SA limits :
lat_max_SA = max(lat_max_list);
lat_min_SA = min(lat_min_list);
lon_max_SA = max(lon_max_list);
lon_min_SA = min(lon_min_list);

%%
% image segmentation per squaere(Kmeans)
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
    im = A; 
    imflat = double(reshape(normalize(im), size(im,1) * size(im,2),3)) ; 
    K = 2; 
    [kIDs, kC] = kmeans(imflat, K, 'Display', 'iter', 'MaxIter', 150, 'Start', 'sample') ;
    imout = reshape(kIDs, size(im,1), size(im,2)) ;
    
    % Make the colors consistent. Assuming that there are less 'lead pixel'
    % in an image compared to 'ice pixel', the image will show darker color
    % for minority pixels : 
    
    if mean(imflat(imout==1)) >= mean(imflat(imout==2))
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
        result_N{i}(:,j) = imout([idx_OLCIunderSAR{i}(:,j)]) ; 
        if length(find(result_N{i}(:,j)==0)) >= ceil(Nclose_pixel/2)
            result_ave{i}(:,j) = 2 ; % LEAD = 2
        elseif length(find(result_N{i}(:,j)==50)) >= ceil(Nclose_pixel/2)
            result_ave{i}(:,j) = 1 ; % ICE = 1
        else 
            result_ave{i}(:,j) = 0 ; % AMBIGUOUS (something is wrong ... )
        end 
    end 
% %     
%     
%     %plot the segmented image on the left 
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
%     if not(exist([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA)]))
%         mkdir(['/Users/ericka/Desktop/Thesis/ValidationImageOutput/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA)])
%     end 
%     saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/image',num2str(i),'_','LON',num2str(lon_min_list(i)),'_',num2str(lon_max_list(i)),'.png'])
%     
%   close all 
   clear Oa08_SA Oa05_SA Oa03_SA Oa15_SA Oa21_SA

end

% we now have SAR_idx_persquare that shows the indicies of SAR of Studyarea
% and results_ave that shows the class of those points

% make the indicies in one array 
for i = 1:N_square 
    add1 = SAR_idx_persquare{i} ; 
    add2 = result_ave{i}.' ; 
    if i ==1 
        SAR_idx_all = [add1] ; 
        result_all = [add2] ; 
    else 
        SAR_idx_all = [SAR_idx_all; add1] ; 
        result_all = [result_all; add2] ; 

    end 
end 

%
% dealing with the overlaps
% making sure that the repeated indicies are the ones that the results will
% be considered. (indicies that appear for the second time)

[SAR_idx_all_nooverlap, ~, ic] = unique(SAR_idx_all, 'last') ; 
result1_all(ic) = result_all ;

save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/result_before_discard.mat'],'result1_all')

    % plot
    figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
    geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
    geoimg.AlphaDataMapping = 'none';
    geoimg.FaceAlpha = 'texturemap';
    alpha(geoimg,double(~isnan(Oa03)))
    hold on 
    scatterm(CS.GEO.LAT(SAR_idx_all(result_all==1)),CS.GEO.LON(SAR_idx_all(result_all==1)), 'filled', 'b'); 
    scatterm(CS.GEO.LAT(SAR_idx_all(result_all==2)),CS.GEO.LON(SAR_idx_all(result_all==2)),  'filled', 'r'); 
    saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/overallresultval1.png'])

    %%
% Removing obvious inaccurate points 
    % now there might be points that are classified into a class that are
    % clearly wrong (refer to the segmentation images)
    % a prompt asking whether there should be a part that should be 'unclassified'

    
discard_lon = [] ;
prompt_YN = 'Are there areas (still) that should be discarded? Y/N : ';
YN = input(prompt_YN, 's') ; 
i = 1; 
while YN =='Y'
   prompt_discardlon = 'Enter range of longitude that should be discarded [min_lon, max_lon] (eg: [172.3, 172.6]): ';
   discard_lon(i,:) = input(prompt_discardlon) ; 
   i = i + 1; 
   YN = input(prompt_YN, 's') ;
end 
 
% find indicies of SAR altimetry that should be discarded 
if ~isempty(discard_lon)
    for i = 1: length(discard_lon(:,1))
        index_to_discard = find(CS.GEO.LON > discard_lon(i,1) & CS.GEO.LON < discard_lon(i,2) ) ;    
        result1_all(ismember(SAR_idx_all_nooverlap ,index_to_discard)==1) = nan ; 
        clear index_to_discard
    end 
end 
% discard (=nan) the results that corresponds to those indicies 
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/result1.mat'],'result1_all')


%% ---------------------------- ValidationScheme 2 -----------------------------------
% ------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------
%looking at sudden changes of radiance along track of altimetry measurements 

%Cropping out irrelevant indicies in the image matrix 
crop_idx_lat = find(lat < lat_min_SA | lat > lat_max_SA ) ; 
crop_idx_lon = find(lon < lon_min_SA | lon > lon_max_SA ) ;
crop_idx = intersect(crop_idx_lat, crop_idx_lon)  ;  % these are the indices that should be cropped !

Oa21(crop_idx) = NaN; 
Oa20(crop_idx) = NaN;
Oa08(crop_idx) = NaN; 
Oa05(crop_idx) = NaN;
Oa03(crop_idx) = NaN;
Oa01(crop_idx) = NaN;

%
%select which band's radiation to  use
rad = Oa03 ;  % Oa08 or the NDSIII
Nclose_pixel = 1; 

[idx_OLCIunderSAR, SAR_SAintersect_idx] = correspondingOLCI(CS, lat, lon, lon_max_SA, lon_min_SA, Nclose_pixel ) ; 
wftosave = CS.SAR.data(:,SAR_SAintersect_idx) ; 
lonSARtosave = CS.GEO.LON(SAR_SAintersect_idx) ; 
latSARtosave = CS.GEO.LAT(SAR_SAintersect_idx) ; 
sigma0 = CS.sigma0(SAR_SAintersect_idx) ; 
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/SAR_idx.mat'],'SAR_SAintersect_idx')
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/sigma0.mat'],'sigma0')
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/wfSA.mat'],'wftosave')
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/lonSAR.mat'],'lonSARtosave')
save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/latSAR.mat'],'latSARtosave')


rad_series = mean(rad(idx_OLCIunderSAR), 1) ; 

% figure
% plot(1:length(rad_series), rad_series)
 max_pks= [];
 max_pks_idx = [];
 min_pks= [];
 min_pks_idx = [];
 [max_pks, max_pks_idx] = findpeaks(rad_series) ; 
 [min_pks, min_pks_idx] = findpeaks(-rad_series) ; 

    % save the image 
    ROOT_imageoutput = '/Users/ericka/Desktop/Thesis/ValidationImageOutput' ; 
    if not(exist([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA)]))
        mkdir(['/Users/ericka/Desktop/Thesis/ValidationImageOutput/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA)])
    end 


figure, findpeaks(rad_series)
saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/radiancechangewithpeaks.png'])


%min peaks are the possible leads (if large enough compared to previous max) 
%for every min peak, the difference to the previous max peak will be
% checked.  
previous_max_idx = [];
previous_max = []; 
next_max = []; 
next_max_idx = [];
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

% threshold is given such that if the dip is larger than that ratio
% compared to the previous max peak, then that index is considered as a
% lead
threshold = 0.015;  
lead_min_pks = []; 
lead_min_pks = find(abs(min_pks./previous_max) < (1- threshold)); 
lead_idx = [];
for i = 1:length(lead_min_pks)
    j = lead_min_pks(i) ; 
    append = previous_max_idx(j)+1:next_max_idx(j)-1 ; % maybe add if statement if the lead is very dark it can be wider.. 
    if ~isnan(append)
        lead_idx = [lead_idx append] ; 
    end 
end 

%%
% Now check the radiance change graph and the peaks found in them 
% if there are very wide leads, there are chances that there will be some 
% maxi/min peaks detected within the dip -> resulting in classifying that
% wide lead as non-lead class. 
% In such cases, by visual inspection, we will specify the indicies that
% should be considered as lead. 

prompt_YN = 'Are there (still) wide leads that are wrongly classified? Y/N : ';
YN = input(prompt_YN, 's') ; 
i = 1; 
while YN =='Y'
   prompt_addlead = 'Enter range of indicies which should be considered as leads [min_idx,max_idx] (eg: [4,0]): ';
   add_lead(i,:) = input(prompt_addlead) ;
   append = [add_lead(i,1) : add_lead(i,2)] ; 
   lead_idx = unique([lead_idx append]) ;
   i = i +1; 
   YN = input(prompt_YN, 's') ;
end 


OLCI_lead_lat = lat([idx_OLCIunderSAR(lead_idx)]) ; 
OLCI_lead_lon = lon([idx_OLCIunderSAR(lead_idx)]) ; 

result2_all = ones(1, length(SAR_SAintersect_idx)) ; 
result2_all(lead_idx) = 2; 

% Removing obvious inaccurate points 
    % now there might be points that are classified into a class that are
    % clearly wrong (refer to the segmentation images)
    % a prompt asking whether there should be a part that should be 'unclassified'

discard_lon2=[];
prompt_YN = 'Are there areas (still) that should be discarded? Y/N : ';
YN = input(prompt_YN, 's') ; 
i = 1; 
while YN =='Y'
   prompt_discardlon = 'Enter range of longitude that should be discarded [min_lon, max_lon] (eg: [172.3, 172.6]): ';
   discard_lon2(i,:) = input(prompt_discardlon) ; 
   i = i +1; 
   YN = input(prompt_YN, 's') ;
end 
 
if ~isempty(discard_lon2)
    % find indicies of SAR altimetry that should be discarded 
    for i = 1: length(discard_lon2(:,1))
        index_to_discard2 = find(CS.GEO.LON > discard_lon2(i,1) & CS.GEO.LON < discard_lon2(i,2) ) ;    
        result2_all(ismember(SAR_idx_all_nooverlap ,index_to_discard2)==1) = nan ; 
        clear index_to_discard2
    end 
end 

% discard (=nan) the results that corresponds to those indicies 

save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/result2.mat'],'result2_all')

   % plot
%     figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
%     geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
%     geoimg.AlphaDataMapping = 'none';
%     geoimg.FaceAlpha = 'texturemap';
%     alpha(geoimg,double(~isnan(Oa03)))
%     hold on 
%     scatterm(CS.GEO.LAT(SAR_idx_all(result2_all==1)),CS.GEO.LON(SAR_idx_all(result2_all==1)), 'filled', 'b'); 
%     scatterm(CS.GEO.LAT(SAR_idx_all(result2_all==2)),CS.GEO.LON(SAR_idx_all(result2_all==2)),  'filled', 'r'); 
%     saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/overallresultval2.png'])

    
% Analyse results obtained by validation 1 and 2: 

lead_confidence = zeros(1,length(result1_all)) ; 
for i = 1: length(result1_all)
    if result1_all(i) == 2 && result2_all(i) == 2 % both show leads
        lead_confidence(i) = 3;
    elseif isnan(result1_all(i)) && result2_all(i) == 2 % unrecognized by val1 but val2 = leads
        lead_confidence(i) = 4;
    elseif isnan(result1_all(i)) && isnan(result2_all(i)) % BOTH NAN
        lead_confidence(i) = nan;           
    elseif result1_all(i) == 1 && result2_all(i) == 1 % both show non-lead
        lead_confidence(i) = 1;
    else 
        lead_confidence(i) = 2;  % one shows lead one shows ice 
    end 
end 

save([ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/lead_confidence.mat'],'lead_confidence')
lead_val_idx = SAR_SAintersect_idx(lead_confidence == 3 | lead_confidence==4 ) ; 
nonlead_val_idx = SAR_SAintersect_idx(lead_confidence == 1 | lead_confidence== 2) ; 


% Analysing waveforms of the result ! 
maxwf = max(CS.SAR.data(:,lead_val_idx)) ; 

length(find(maxwf < 1000)) 
length(find(maxwf < 5000)) 

figure 
subplot(2,3,1)
plot(CS.SAR.data(:,lead_val_idx(maxwf < 1000)))
subplot(2,3,2)
plot(CS.SAR.data(:,lead_val_idx(maxwf > 1000 & maxwf < 5000)))
subplot(2,3,3)
plot(CS.SAR.data(:,lead_val_idx(maxwf > 5000 & maxwf < 10000)))
subplot(2,3,4)
plot(CS.SAR.data(:,lead_val_idx(maxwf > 10000 & maxwf < 50000)))
subplot(2,3,5)
plot(CS.SAR.data(:,lead_val_idx(maxwf > 50000 & maxwf < 100000)))
subplot(2,3,6)
plot(CS.SAR.data(:,lead_val_idx(maxwf > 100000)))
saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/final_leadwaveforms.png'])


figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
geoimg.AlphaDataMapping = 'none';
geoimg.FaceAlpha = 'texturemap';
alpha(geoimg,double(~isnan(Oa03)))
hold on 
scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==1)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==1)), 'filled', 'b') ;
scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==2)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==2)), 'filled', 'b') ;
scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==3)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==3)), 'filled', 'r') ;
scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==4)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==4)), 'filled', 'r') ;
saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/finalmapresult.png'])

