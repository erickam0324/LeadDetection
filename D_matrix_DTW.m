function D = D_matrix_DTW(waveform)
    % function to calculate the D matrix based on Dynanic Time Warping (DTW)
    % input is a (128 x n) waveform data  

    nwf = length(waveform(1,:)) ; 
    
    %% Assemble
    D= zeros(nwf,nwf);
    for i = 1 : nwf 
        for j =  i+1:nwf
        D(i,j) = dtw(waveform(:,i), waveform(:,j)) ; 
        end 
    end 
    D = D + D.'; 
    
end
