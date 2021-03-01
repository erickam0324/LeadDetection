%% 
function labels = fun_DBSCAN_classification(waveform, epsilon, Minpts )
global D 
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% X here is the waveform !!! 

[class, isnoise, D] = DBSCAN(waveform,epsilon,Minpts); 
nclass = length(unique(class)); % number of class that was found by DBSCAN

class = class+1; 

figure 
for i = 1:nclass
     subplot(ceil(nclass/5),5,i)
     ind = class==i; 
     plot(waveform(:,ind))
     hold on
end
    sgtitle("Multi-parameter DBSCAN cluster")


    %% Asks input from the user 
    prompt_lead = 'Which graph number can be classified as leads? If there are multiple, enter as an array (eg; [3,4,5]): ';
    lead_k = input(prompt_lead) ; 
%     prompt_ice = 'Which graph number can be classified as ice? If there are multiple, enter as an array (eg; [3,4,5]): ';
%     ice_k = input(prompt_ice) ; 
%     prompt_ocean = 'Which graph number can be classified as ocean? If there are multiple, enter as an array (eg; [3,4,5]): ';
% %     ocean_k = input(prompt_ocean) ; 
%     ambiguous_k = setdiff(1:k_final, lead_k) ; 
%     ambiguous_k = setdiff(ambiguous_k, ice_k) ; 
%     ambiguous_k = setdiff(ambiguous_k, ocean_k) ; 
    ice_k = setdiff(1:nclass, lead_k) ; 
    ocean_k = [] ; 
    ambiguous_k = [] ; 
    
    labels = nan(1,numel(waveform(1,:))) ; %ice = 3
    leadindex = ismember(class,lead_k)==1 ; 
    iceindex = ismember(class,ice_k)==1; 
    labels(leadindex) = 3;
    labels(iceindex) = 1; 
    
    

    function [IDX, isnoise, D]=DBSCAN(X,epsilon,MinPts)
        C=0;

        n=size(X,2);
        IDX=zeros(n,1);

        %D=load('D_DTW.mat') ;
        %D = D.D ; 
        D = D_matrix(X); 

        visited=false(n,1);
        isnoise=false(n,1);

        for i=1:n
            if ~visited(i)
                visited(i)=true;

                Neighbors=RegionQuery(i);
                if numel(Neighbors)< MinPts
                    % X(i,:) is NOISE
                    isnoise(i)=true;
                else
                    C=C+1;
                    ExpandCluster(i,Neighbors,C);
                end

            end

        end

        function ExpandCluster(i,Neighbors,C)
            IDX(i)=C;

            k = 1;
            while true
                j = Neighbors(k);

                if ~visited(j)
                    visited(j)=true;
                    Neighbors2=RegionQuery(j);
                    if numel(Neighbors2)>=MinPts
                        Neighbors=[Neighbors Neighbors2];   %#ok
                    end
                end
                if IDX(j)==0
                    IDX(j)=C;
                end

                k = k + 1;
                if k > numel(Neighbors)
                    break;
                end
            end
        end

        function Neighbors=RegionQuery(i)
            Neighbors=find(D(i,:)<=epsilon);
        end

    end


end

