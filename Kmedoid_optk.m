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
%wf   = NORMfactor .* CS.SAR.data;
wf = CS.SAR.data;

%newwf = [wf(:,1) wf(:,2) wf(:,3) wf(:,4) wf(:,5) wf(:,6) wf(:,4998) wf(:,4999) wf(:,5000) wf(:,5001)] ; 
newwf = wf(:,1:3000) ; 
nwf = length(newwf(1,:)) ; 

for i = 1:nwf
    for j = 1: nwf 
        D(i,j) = sum((newwf(:,i) - newwf(:,j)).^2); 
    end 
end 

% if k is too large for selected nwf, the index overlaps and in the end we will end up with less clusters 
n= nwf;

itr=0; 
k_list = 2:2:50 ; 
for k = k_list
    itr = itr+1; 
    label = ceil(k*rand(1,n));
    last = zeros(1,n);
    while any(label ~= last)                            % continue until convergence 
        [~,~,last(:)] = unique(label);                  % remove empty clusters
        [~, index] = min(D*sparse(1:n,last,1),[],1);    % find k medoids
        [val, label] = min(D(index,:),[],1);            % assign labels
    end
    energy(itr) = sum(val);
    k_final(itr) = length(index); 
end 

figure
plot(k_list, energy)
xlabel('Number of Kmedoid cluster')
ylabel('Energy (total distances to the medoid in all clusters)')
title("Cluster number vs performance of Kmedoid clustering with n=" + nwf + " waveform data")

