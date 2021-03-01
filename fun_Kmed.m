%% FUNCTION FOR KMED + KNN CLASSIFICATION 
% INPUT : waveform (n x m) matrix 
% OUTPUT : (m) vector w/ classification identifier number per waveform 

function [energy, k_final] = fun_Kmed(wfpool, k)
    %% Read S3 L1b data and compute normalized waveforms 

    nwfpool = length(wfpool(1,:)) ; 

    [~, outlier_all, idx_all_nooutliers, wf_all_nooutliers] = D_matrix(wfpool); 

    [Dall_nooutliers, ~] = D_matrix(wf_all_nooutliers); 
    nwf_all_nooutliers = size(Dall_nooutliers,1) ; 

    %% Kmedoid (setting up training data)

    % Select fraction from the Kmedoid clusters to create a KNN training set 
    n_train =  ceil(0.3*nwf_all_nooutliers); 
    n_test = nwf_all_nooutliers - n_train  ; 
    train_idx = sort(randperm(nwf_all_nooutliers,n_train)) ; % selecting indicies from non outlier waveform 
    test_idx = setdiff(1:nwf_all_nooutliers, train_idx) ; 
    train_idx_all = idx_all_nooutliers(train_idx) ; % selected wavforms in terms of the very original index (wfpool)
    test_idx_all = setdiff(idx_all_nooutliers,train_idx_all) ; 

    wf_train = wfpool(:,train_idx_all) ; 
    wf_test =  wfpool(:,test_idx_all) ; 


    [D_km, ~] = D_matrix(wf_train); 
    
    label_train = ceil(k*rand(1,n_train));

    last = zeros(1,n_train);
    while any(label_train ~= last)                            % continue until convergence 
        [~,~,last(:)] = unique(label_train);                  % remove empty clusters
        [~, index] = min(D_km*sparse(1:n_train,last,1),[],1);    % find k medoids
        [val, label_train] = min(D_km(index,:),[],1);              % assign labels
    end
    energy = sum(val);
    k_final = length(index) ; 

%     figure 
%     for i = 1:k_final
%         subplot(ceil(k_final/5),5,i)
%         for j = 1:length(find(label_train==i))
%             ind = find(label_train==i); 
%             plot(wf_train(:,ind(j)))
%             hold on
%         end     
%     end
% 
%     sgtitle("Multi-parameter Kmedoid cluster with n=" + nwf_all_nooutliers + " (after oultlier rejection) waveform data")


end 

 


