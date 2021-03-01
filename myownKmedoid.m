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

for i = 1:nwf
    for j = 1: nwf 
        D(i,j) = sum((newwf(:,i) - newwf(:,j)).^2); 
    end 
end 

k = 20;  % if k is too large for selected nwf, the index overlaps and in the end we will end up with less clusters 
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
    subplot(4,3,i)
    for j = 1:length(find(label==i))
        ind = find(label==i); 
        plot(newwf(:,ind(j)))
        hold on
    end     
end

sgtitle("Euclidean distance Kmedoid cluster with n=" + nwf + " waveform data")

