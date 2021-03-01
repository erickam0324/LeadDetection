%% FUNCTION FOR KMED + KNN CLASSIFICATION 
% INPUT : waveform (n x m) matrix 
% OUTPUT : (m) vector w/ classification identifier number per waveform 

function LABEL_FINAL_01234 = fun_Kmed_KNN_classification(wfpool, k)
    %% Read S3 L1b data and compute normalized waveforms 

    nwfpool = length(wfpool(1,:)) ; 

    [~, outlier_all, idx_all_nooutliers, wf_all_nooutliers] = D_matrix(wfpool); 

    
    [Dall_nooutliers, ~] = D_matrix(wf_all_nooutliers); 
    nwf_all_nooutliers = size(Dall_nooutliers,1) ; 

    %% Kmedoid (setting up training data)
    k_knn = 5 ; 

    % Select fraction from the Kmedoid clusters to create a KNN training set 
    n_train =  ceil(0.3*nwf_all_nooutliers); 
    n_test = nwf_all_nooutliers - n_train  ; 
    train_idx = sort(randperm(nwf_all_nooutliers,n_train)) ; % selecting indicies from non outlier waveform 
    test_idx = setdiff(1:nwf_all_nooutliers, train_idx) ; 
    train_idx_all = idx_all_nooutliers(train_idx) ; % selected wavforms in terms of the very original index (wfpool)
    test_idx_all = setdiff(idx_all_nooutliers,train_idx_all) ; 

    wf_train = wfpool(:,train_idx_all) ; 
    wf_test =  wfpool(:,test_idx_all) ; 


    [D_km] = D_matrix(wf_train); 


    % if k is too large for selected nwf, the index overlaps and in the end we will end up with less clusters 

    label_train = ceil(k*rand(1,n_train));

    last = zeros(1,n_train);
    while any(label_train ~= last)                            % continue until convergence 
        [~,~,last(:)] = unique(label_train);                  % remove empty clusters
        [~, index] = min(D_km*sparse(1:n_train,last,1),[],1);    % find k medoids
        [val, label_train] = min(D_km(index,:),[],1);              % assign labels
    end
    energy = sum(val);
    k_final = length(index) 

    figure 
    for i = 1:k_final
        subplot(ceil(k_final/5),5,i)
        for j = 1:length(find(label_train==i))
            ind = find(label_train==i); 
            plot(wf_train(:,ind(j)))
            hold on
        end     
    end

    sgtitle("Multi-parameter Kmedoid cluster with n=" + nwf_all_nooutliers + " (after oultlier rejection) waveform data")


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

    ice_k = setdiff(1:k_final, lead_k) ; 
    ocean_k = [] ; 
    ambiguous_k = [] ; 
    
    labels = label_train ; % label of training data
    %labels_test = label(test_idx) ; % label of test data
    labels_test =zeros(1,n_test) ;

    %% KNN part 

    % initialize matrices
    ed=zeros(n_test, n_train); %ed: (MxN) euclidean distances 
    ed_ind=zeros(n_test, n_train);
    k_nn=zeros(n_test,k_knn); %k-nearest neighbors for testing sample (Mxk)
    predicted_labels=zeros(n_test,1);

    %----------------------------------------------------------------
    %-------MAKE THIS PART MORE EFFICIENT !!!!!!! -------------------
    %----------------------------------------------------------------
    % % identify and get rid of outliers also for the testing data. 
    % 
    % [~, outlier_test, idxtest_nooutliers, wftest_nooutliers] = D_matrix(wf_test); 
    % n_test_nooutliers = n_test - length(outlier_test) ; 
    % outlier_test_all = test_idx(outlier_test) ; 

    % Make the distance matrix by reorganizing the original distance matrix D
    for i = 1 :n_test
        current_test_idx = test_idx(i) ; 
        for j = 1 : n_train
        current_train_idx = train_idx(j) ;
        ed(i,j) = Dall_nooutliers(current_test_idx, current_train_idx) ; 
        end
        [ed(i,:),ed_ind(i,:)]=sort(ed(i,:)) ;
    end

    %find the nearest k for each data point of the testing data
    k_nn=ed_ind(:,1:k_knn);
    nn_index=k_nn(:,1);
    %get the majority vote 
    for i=1:size(k_nn,1)
        options=unique(labels(k_nn(i,:)'));
        max_count=0;
        max_label=0;
        for j=1:length(options)
            L=length(find(labels(k_nn(i,:)')==options(j)));
            if L>max_count
                max_label=options(j);
                max_count=L;
            end
        end
        predicted_labels(i)=max_label;
    end

    label_final = nan(1,nwfpool) ; 
    label_final(train_idx_all) = labels; 
    label_final(test_idx_all) = predicted_labels; %sorted back to the original order of SAR waveform 

    leadindex = ismember(label_final,lead_k)==1 ; 
    iceindex = ismember(label_final,ice_k)==1; 
    oceanindex = ismember(label_final,ocean_k)==1; 
    ambiguousindex = ismember(label_final,ambiguous_k)==1;

    % AMBIGUOUS = 0
    % ICE = 1
    % LEAD = 3
    % OCEAN = 2
    % OUTLIERS = 4

    LABEL_FINAL_01234 = nan(1, nwfpool); 
    LABEL_FINAL_01234(ambiguousindex) = 0;
    LABEL_FINAL_01234(iceindex) = 1;
    LABEL_FINAL_01234(leadindex) = 3; 
    LABEL_FINAL_01234(oceanindex) = 2;
    LABEL_FINAL_01234(outlier_all) = 4;
end 

 


