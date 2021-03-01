%% Function to find corresponding latitude for given longitude inputs 

function [lat_min_SA, lat_max_SA] = correspondingLAT(CS, lon_min_SA, lon_max_SA, margin) 
    lat_max_SA = CS.GEO.LAT(find(CS.GEO.LON < lon_max_SA, 1 )) + margin ; 
    lat_min_SA = CS.GEO.LAT(find(CS.GEO.LON > lon_min_SA, 1, 'last' )) - margin;
end 