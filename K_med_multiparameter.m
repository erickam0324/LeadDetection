%% multiparameter k medoid
% First we need waveform data! 
% Select 10 different waveform from sen

clear all; 
LoadCommonSettings_ericka;

rng(19970324)

SAT = 'S3A' ; 
defval('DOM',[60 90; -180 180])               %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                 %Retracker to be used
defval('SetRetr',{})
defval('DEM',[])    
defval('IDXpoi',[])          

FName = 'S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';

%% Read S3 L1b data and compute normalized waveforms 

[~,CS] = S3_L1b_read_ericka(FName)   ;
DDAcf   = DDA_ConfigFile(SAT,'SAR');
Wf_Norm_Aw = 1;      
%NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
%wf   = NORMfactor .* CS.SAR.        c   v     v cv c c c c data;
wf = CS.SAR.data;

%newwf = [wf(:,1) wf(:,2) wf(:,3) wf(:,4) wf(:,5) wf(:,6) wf(:,4998) wf(:,4999) wf(:,5000) wf(:,5001)] ; 
newwf = wf(:,3000:6000) ; 
nwf = length(newwf(1,:)) ; 

%% Distances 
for i = 1:nwf
    for j = 1: nwf 
        D_d(i,j) = sum((newwf(:,i) - newwf(:,j)).^2); 
    end 
end 
D_d = D_d./ max(D_d) ; 

%% MAX 
[MAX_value, MAX_bin] = max(newwf,[],1) ;

for i = 1: nwf 
    for j = 1:nwf
    D_MAX(i,j) = (MAX_value(i) - MAX_value(j))^2; 
    end 
end 
D_MAX = D_MAX./ max(D_MAX) ; 

%% Waveform width 
for i = 1:nwf 
    ww(i) = length(find (newwf(:,i)> 0.01*max(newwf(:,i)))) ; 
end
for i = 1: nwf 
    for j = 1:nwf
        D_ww(i,j) = (ww(i) - ww(j))^2; 
    end 
end 
D_ww = D_ww./ max(D_ww) ; 

%% LeS 
for i = 1:nwf 
    MAX_12_bin_list = find(newwf(:,i)> 0.125*MAX_value(i)) ; 
    LeS(i) = MAX_bin(i) - MAX_12_bin_list(1) ; 
end
for i = 1: nwf 
    for j = 1:nwf
    D_LeS(i,j) = (LeS(i) - LeS(j))^2; 
    end 
end 
D_LeS = D_LeS./ max(D_LeS) ; 

%% TeS 
for i = 1:nwf 
    MIN_12_bin_list = find(newwf(:,i)< 0.125*MAX_value(i)) ; 
    TeS(i) = MAX_bin(i) - min(MIN_12_bin_list(MIN_12_bin_list > MAX_bin(i))); 
end
for i = 1: nwf 
    for j = 1:nwf
    D_TeS(i,j) = (TeS(i) - TeS(j))^2; 
    end 
end 
D_TeS = D_TeS./ max(D_TeS) ;

%% 
D = D_TeS + D_LeS + D_ww + D_MAX; 
k = 30;  % if k is too large for selected nwf, the index overlaps and in the end we will end up with less clusters 
n= nwf;
label = ceil(k*rand(1,n));

last = zeros(1,n);
while any(label ~= last)                            % continue until convergence 
    [~,~,last(:)] = unique(label);                  % remove empty clusters
    [~, index] = min(D*sparse(1:n,last,1),[],1);    % find k medoids
    [val, label] = min(D(index,:),[],1);            % assign labels
end
energy = sum(val);

k_final = length(index) ; 

figure 
for i = 1:k_final
    subplot(3,4,i)
    for j = 1:length(find(label==i))
        ind = find(label==i); 
        plot(newwf(:,ind(j)))
        hold on
    end     
end
sgtitle("Multi-parameter Kmedoid cluster with n=" + nwf + " waveform data")


