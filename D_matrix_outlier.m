
%INPUT : 
% waveform  = waveform data 
% outlier_yesORno = if "yes" then the function will look for outliers and
% compute D matrix without outlier indices 


function [D, outlier_index, idx_nooutliers, wf_nooutliers] = D_matrix(waveform)
    nwf = length(waveform(1,:)) ; 
    %% Identifying outlier waveforms 
    outlier = zeros(1,nwf); 
    %% MAX 
    [MAX_value, MAX_bin] = max(waveform,[],1) ;
    
    %% Waveform width 
    for i = 1:nwf 
        ww(i) = length(find (waveform(:,i)> 0.01*max(waveform(:,i)))) ;
    end
    
    %% LeS 
    for i = 1:nwf 
        MAX_12_bin_list = find(waveform(:,i)> 0.125*MAX_value(i)) ; 
        if MAX_bin(i) == 1
            outlier(i) = 1; 
        elseif isempty(min(MAX_12_bin_list(MAX_12_bin_list < MAX_bin(i))))
            if  waveform(MAX_bin(i)-1,i) < 0.125*MAX_value(i)
                outlier(i) = 0;  
                LeS(i) = 0; 
            else
                outlier(i) = 1; 
                LeS(i) = 0; 
            end 
        else
            LeS(i) = MAX_bin(i) - MAX_12_bin_list(1) ; 
        end 
    end
    
    
    %% TeS
    TeS = zeros(1,nwf);
    for i = 1:nwf 
        MIN_12_bin_list = find(waveform(:,i)< 0.125*MAX_value(i)) ; 
        if isempty(MIN_12_bin_list) 
            TeS(i) = 129          ;         % some outlier waveforms have very shallow slope that never reaches 12.5% -> give very shallow value for TeS
            outlier(i) = 1; 
        elseif  isempty(min(MIN_12_bin_list(MIN_12_bin_list > MAX_bin(i))))
            TeS(i) = mean(TeS)         ;    % look at i = 8282, the waveform is incomplete 
            outlier(i) = 1; 
        else 
            TeS(i) = min(MIN_12_bin_list(MIN_12_bin_list > MAX_bin(i))) - MAX_bin(i); 
        end
    end 
    
    
    %% Dealing with outliers 
    outlier_index = find(outlier == 1) ; % find indicies with outliers 
    idx_nooutliers = setdiff([1:nwf], outlier_index) ; % remove outlier indicies 
    wf_nooutliers = waveform(:,idx_nooutliers) ; 
   
        %% Assembling the D_matrix without outliers 
        if ~isempty(outlier_index)
            MAX_value(outlier_index) = []; 
            ww(outlier_index) = [];     
            LeS(outlier_index) = []; 
            TeS(outlier_index) = []; 
        end 
        new_nwf = length(wf_nooutliers(1,:)) ; 
        
    D_MAX = zeros(new_nwf,new_nwf);
    for i = 1 : new_nwf 
        for j = 1+i:new_nwf
        D_MAX(i,j) = (MAX_value(i) - MAX_value(j))^2; 
        end 
    end 
    D_MAX = D_MAX + D_MAX.';
    D_MAX = D_MAX./ max(D_MAX) ; 
    D_MAX = D_MAX.*2 ; 
    disp('DMAX done')
    
    D_ww = zeros(new_nwf,new_nwf);
    for i = 1: new_nwf
        for j = 1+i:new_nwf
            D_ww(i,j) = (ww(i) - ww(j))^2; 
        end 
    end 
    D_ww = D_ww + D_ww.' ; 
    D_ww = D_ww./ max(D_ww) ; 
    disp('Dww done')

    D_LeS = zeros(new_nwf,new_nwf) ;
    for i = 1: new_nwf 
        for j = i+1:new_nwf
        D_LeS(i,j) = (LeS(i) - LeS(j))^2; 
        end 
    end 
    D_LeS = D_LeS + D_LeS.'; 
    D_LeS = D_LeS./ max(D_LeS) ; 
    disp('Dles done')
    
    D_TeS = zeros(new_nwf,new_nwf) ; 
    for i = 1: new_nwf
        for j = i+1:new_nwf
        D_TeS(i,j) = (TeS(i) - TeS(j))^2; 
        end 
    end 
    D_TeS = D_TeS + D_TeS.'; 
    D_TeS = D_TeS./ max(D_TeS) ;
    disp('Dtes done')
    
    %% 
    D = D_TeS + D_LeS + D_ww + D_MAX; 
   
end
