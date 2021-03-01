%% READ S3 L1b data
clear all; 
LoadCommonSettings_ericka;


SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          

for ii = 2
    OLCIpath = '/Users/ericka/Desktop/Thesis/OLCIdata';
    S1 = dir(fullfile(OLCIpath,'*'));
    NamesOLCIdata = setdiff({S1([S1.isdir]).name},{'.','..'}); % list of subfolders of D.

    SARpath = '/Users/ericka/Desktop/Thesis/SARdata';
    S2 = dir(fullfile(SARpath,'*'));
    NamesSARdata = setdiff({S2([S2.isdir]).name},{'.','..'}); % list of subfolders of D.

    FName = NamesSARdata{ii} ; 
    FName = fullfile(SARpath, FName) ;
    [~,CS, NORMfactor] = S3_L1b_read_ericka(FName)   ;



    %% Code to raed OLCI GEO TIFF IMAGE 

    % Locate OLCI Geo Tiff image : 
    SEN3file = NamesOLCIdata{ii} ; 
    date = SEN3file(17:24) ; 
    SEN3file = fullfile(OLCIpath, SEN3file) ;


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
    if ~isempty(find(lon>0 & lon<1, 1)) && ~isempty(find(lon<0 & lon> -1, 1)) && ~isempty(find(lon<180 & lon> 179, 1)) && ~isempty(find(lon<-179 & lon> -180, 1))
        % this is the rare case when olci image spans for the whole
        % longitde range (has values between -1 and 1 but also 179 & -179)
        % use 0 to 360 coordinate system (change OLCI lon values)
        coordinate = '0to360';
            for i = 1:length(lon(:,1))
                for j = 1:length(lon(1,:))
                    if lon(i,j) < 0 
                        lon(i,j) = lon(i,j) + 360;
                    end 
                end
            end        
    elseif ~isempty(find(lon>0 & lon<1, 1)) && ~isempty(find(lon<0 & lon> -1, 1))
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

end 
