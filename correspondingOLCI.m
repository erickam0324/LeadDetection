%% Function for finding corresponding pixel to the given SAR altimetry track 
% Input : SAR GEO data, OLCI GEO data, study area (min/max of lon/lat), #of closest
% pixels
% Output : OCLI idx & OLCI values 


function [output, SAR_SAintersect_idx] = correspondingOLCI(CS, lat, lon, lon_max_SA, lon_min_SA, Nclose_pixel )

% SAR's corresponding latitude given the Study area longitude 
[lat_min_SA, lat_max_SA] = correspondingLAT(CS, lon_min_SA, lon_max_SA, 0.05) ; 

% Find index of SAR data which corresponds to the Study area given above 
LONindex_SA = find(CS.GEO.LON > lon_min_SA & CS.GEO.LON < lon_max_SA ) ;
LATindex_SA =  find(CS.GEO.LAT > lat_min_SA & CS.GEO.LAT <lat_max_SA ) ; 
SAR_SAintersect_idx = intersect(LONindex_SA, LATindex_SA) ;  % indicies of SAR altimetry in the SA
LON_SAR_SA = CS.GEO.LON(SAR_SAintersect_idx); % lon values of SAR altimetry in the SA
LAT_SAR_SA = CS.GEO.LAT(SAR_SAintersect_idx); % lat values of SAR altimetry in the SA

% initialize 
idx_OLCIunderSAR = [];
    for i = 1:length(SAR_SAintersect_idx)
        perturb = 0.03 ; 
        % look for olci idx that are +- pertubation of both lon and lat of the SAR point 
        olci_search_lonidx = find(lon > LON_SAR_SA(i)-perturb & lon < LON_SAR_SA(i)+perturb ) ;
        olci_search_latidx =  find(lat > LAT_SAR_SA(i)-perturb & lat < LAT_SAR_SA(i)+perturb ) ;
        olci_searchintersect_idx = intersect(olci_search_lonidx, olci_search_latidx) ; 

        % compute distnace between these selcted points to the SAR altimetry point
        % of interest : 
        distance = abs(lon(olci_searchintersect_idx) - LON_SAR_SA(i)) + abs(lat(olci_searchintersect_idx) - LAT_SAR_SA(i)) ; 
        [~,I] = mink(distance,Nclose_pixel); 
        
        idx_OLCIunderSAR(:,i) = olci_searchintersect_idx(I) ; % indicies of olci that are 'under' SAR. 
    end 
    output = idx_OLCIunderSAR ; 
end 
