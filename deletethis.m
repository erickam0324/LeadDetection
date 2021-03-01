%% Delete this soon 

clear all; 
LoadCommonSettings_ericka;
rng(19970324)

D = '/Users/ericka/Desktop/Thesis/ValidationImageOutput';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for k = 1:numel(N)
    SAR_idx{k} = load(fullfile(D,N{k},'SAR_idx.mat'));
end

for ii = 1: 2

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
   

    ROOT_imageoutput = '/Users/ericka/Desktop/Thesis/ValidationImageOutput' ; 

    sigma0 = CS.sigma0(SAR_idx{ii}.SAR_SAintersect_idx) ; 
    save([ROOT_imageoutput,'/',N{ii},'/sigma0.mat'],'sigma0')
 end


%%
% GEOfile = fullfile(SEN3file, "geo_coordinates.nc") ; 
% Oa01 = ncread(fullfile(SEN3file, "Oa01_radiance.nc"), "Oa01_radiance");
% Oa03 = ncread(fullfile(SEN3file, "Oa03_radiance.nc"), "Oa03_radiance");
% Oa05 = ncread(fullfile(SEN3file, "Oa05_radiance.nc"), "Oa05_radiance");
% Oa08 = ncread(fullfile(SEN3file, "Oa08_radiance.nc"), "Oa08_radiance");
% Oa15 = ncread(fullfile(SEN3file, "Oa15_radiance.nc"), "Oa15_radiance");
% Oa21 = ncread(fullfile(SEN3file, "Oa21_radiance.nc"), "Oa21_radiance");
% lat = ncread(GEOfile, "latitude");
% lon = ncread(GEOfile, "longitude");
% 
% % LOOKING AT WHERE OLCI IMAGE IS LOCATED. 
% % if olci is has coordinates that are crossing values at 0deg, the
% % coordinates should be expressed in -180 to 180. Otherwise, 0 to 360 
% if ~isempty(find(lon>0 & lon<1, 1)) && ~isempty(find(lon<0 & lon> -1, 1))
%     % use -180 to 180 coordinate system (change SAR lon values)
%     coordinate = '-180to180' ; 
%         for i = 1:length(CS.GEO.LON)
%             if CS.GEO.LON(i) >180
%                CS.GEO.LON(i) =  CS.GEO.LON(i) - 360;
%             end 
%         end
% else 
%     % use 0 to 360 coordinate system (change OLCI lon values)
%     coordinate = '0to360';
%         for i = 1:length(lon(:,1))
%             for j = 1:length(lon(1,:))
%                 if lon(i,j) < 0 
%                     lon(i,j) = lon(i,j) + 360;
%                 end 
%             end
%         end
% end 
% 
% A_col(:,:,1) = Oa08;  % red
% A_col(:,:,2) = Oa05;  % green
% A_col(:,:,3) = Oa03;  % blue
% 
% % Load SAR altimetery files 
% SAT = 'S3A' ; 
% defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
% defval('Retracker','SAMOSA2')                 %Retracker to be used
% defval('SetRetr',{})
% defval('DEM',[])    
% defval('IDXpoi',[])          
% 
% 
% % Find index of SAR data which corresponds to the image coordinates 
% LONindex = find(CS.GEO.LON > min(lon(:)) & CS.GEO.LON < max(lon(:))) ;
% LATindex =  find(CS.GEO.LAT > min(lat(:)) & CS.GEO.LAT < max(lat(:)) ); 
% OLCI_intersec_idx = intersect(LONindex, LATindex) ; 
% LON_SAR = CS.GEO.LON(OLCI_intersec_idx);
% LAT_SAR = CS.GEO.LAT(OLCI_intersec_idx);
% 


% %% ---------------------------- Study area selection ---------------------------------
% % ------------------------------------------------------------------------------------
% % ------------------------------------------------------------------------------------
% % 
% % enter following lon MAX MIN showing the study area. These values may change
% % slightly depending on the sizes of the small squares that will be used 
% % for image segmentation analysis. 
% lon_min_SA = -3.4; 
% lon_max_SA = 14;
% lat_min_SA = min(CS.GEO.LAT(find(CS.GEO.LON > lon_min_SA & CS.GEO.LON < lon_max_SA )))-0.05 ; 
% lat_max_SA = max(CS.GEO.LAT(find(CS.GEO.LON > lon_min_SA & CS.GEO.LON < lon_max_SA )))+0.05 ; 
% 
% 
% % count how many squares we need to analyse 
% lon_square = 0.6 ; % size (degree of longitude) of a single square
% patch_size = 0.1 ; % this is the size (degree of longitude) that two squares overlap such that validation scheme can be carried out for all points (if 0. the edge points might not have enough data point for validation)
% N_square = ceil((lon_max_SA - lon_min_SA) / (lon_square-patch_size)) ; 
% 
% % Make list of max/min lon/lat per square with overlaps 
% for i = 1:N_square
%    lon_max_list(i) = lon_max_SA - (lon_square-patch_size)*(i-1); 
%    lon_min_list(i) = lon_max_list(i) - lon_square ; 
%    lat_max_list(i) = max(CS.GEO.LAT(find(CS.GEO.LON > lon_min_list(i) & CS.GEO.LON < lon_max_list(i) ))) +0.05 ; 
%    lat_min_list(i) = min(CS.GEO.LAT(find(CS.GEO.LON > lon_min_list(i) & CS.GEO.LON < lon_max_list(i) ))) -0.05 ; 
% 
% end 
% 
% % New SA limits :
% lat_max_SA = max(lat_max_list);
% lat_min_SA = min(lat_min_list);
% lon_max_SA = max(lon_max_list);
% lon_min_SA = min(lon_min_list);





% %%
% D = '/Users/ericka/Desktop/Thesis/ValidationImageOutput';
% S = dir(fullfile(D,'*'));
% N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
% for ii = 1:numel(N)
%     LeadConfidence{ii} = load(fullfile(D,N{ii},'lead_confidence.mat')); % improve by specifying the file extension. 
%     wfSA{ii} = load(fullfile(D,N{ii},'wfSA.mat'));
%     sigma0{ii} = load(fullfile(D,N{ii},'sigma0.mat'));
% end
% %
% % put all data together 
% wf_all = []; 
% conf_all = [];
% sigma0_all = [];
% ii = 7 ; 
% wf_all = [wf_all wfSA{1,ii}.wftosave]; 
% sigma0_all = [sigma0_all sigma0{1,ii}.sigma0tosave.']; 
% conf_all = [conf_all LeadConfidence{1,ii}.lead_confidence] ; 
% 
%     
% lead_confidence = conf_all ; 
%     
% figure, worldmap([lat_min_SA lat_max_SA], [lon_min_SA lon_max_SA])
% geoimg=geoshow(lat, lon, A_col/256,'DisplayType','image');  
% geoimg.AlphaDataMapping = 'none';
% geoimg.FaceAlpha = 'texturemap';
% alpha(geoimg,double(~isnan(Oa03)))
% hold on 
% scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==1)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==1)), 'filled', 'b') ;
% scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==2)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==2)), 'filled', 'y') ;
% scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==3)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==3)), 'filled', 'r') ;
% scatterm(CS.GEO.LAT(SAR_SAintersect_idx(lead_confidence==4)),CS.GEO.LON(SAR_SAintersect_idx(lead_confidence==4)), 'filled', 'r') ;
% saveas(gcf, [ROOT_imageoutput,'/',date,'_LON',num2str(lon_min_SA),'_',num2str(lon_max_SA),'/finalmapresult.png'])
% 
% %%
% 
% lead_val_idx = SAR_SAintersect_idx(lead_confidence == 3 | lead_confidence==4 ) ; 
% nonlead_val_idx = SAR_SAintersect_idx(lead_confidence == 1 | lead_confidence== 2) ; 
% 
% 
% % Analysing waveforms of the result ! 
% maxwf = max(CS.SAR.data(:,lead_val_idx)) ; 
% 
% 
% figure 
% subplot(2,3,1)
% plot(CS.SAR.data(:,lead_val_idx(maxwf < 1000)))
% subplot(2,3,2)
% plot(CS.SAR.data(:,lead_val_idx(maxwf > 1000 & maxwf < 5000)))
% subplot(2,3,3)
% plot(CS.SAR.data(:,lead_val_idx(maxwf > 5000 & maxwf < 10000)))
% subplot(2,3,4)
% plot(CS.SAR.data(:,lead_val_idx(maxwf > 10000 & maxwf < 50000)))
% subplot(2,3,5)
% plot(CS.SAR.data(:,lead_val_idx(maxwf > 50000 & maxwf < 100000)))
% subplot(2,3,6)
% plot(CS.SAR.data(:,lead_val_idx(maxwf > 100000)))
