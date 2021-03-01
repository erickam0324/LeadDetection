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
wfpool = CS.SAR.data;
nwfpool = length(wfpool(1,:)) ; 

%% Selection of wavwforms
% Waveforms are selected for the given waveform numbers 
% first index is selected randomly, and from there, a consequetive number
% of waveforms are selected. 

set_nwf = [1000,2000,3000,4000,5000,6000] ; 
for i = 1:length(set_nwf) 
    first_wf(i) = randi([1 nwfpool-set_nwf(i)]) ; 
    newwf{i} = wfpool(:, (first_wf(i):(first_wf(i)+set_nwf(i)-1))) ; 
end 


k_list = 1:2:50 ; 
energy_list = nan(length(set_nwf), length(k_list)) ; 
for i = 1:length(set_nwf) 
    for j = 1:length(k_list)
        k = k_list(j) ; 
        [energy, k_final] = fun_Kmed(newwf{i}, k) ; 
        energy_list(i,j) = energy ; 
    end 
end 

for i = 1:length(set_nwf) 
[elbowpoint, ep_idx] = knee_pt(energy_list(i,:), k_list) ; 
ep_list(i) = elbowpoint ; 
epidx_list(i) = ep_idx ; 
end 

figure
for i = 1:length(set_nwf) 
    txt = ['n = ',num2str(set_nwf(i))];
    plot(k_list, energy_list(i,:),'DisplayName',txt)
    hold on 
end 
hold on 

for i = 1:length(set_nwf) 
    txt2 = ['e.p. at k = ',num2str(ep_list(i)), ', n = ', num2str(set_nwf(i))];
    plot(ep_list(i), energy_list(i,epidx_list(i)),'*','DisplayName',txt2)
end 
legend show


% figure
% plot(k_list, energy)
% xlabel('Number of Kmedoid cluster')
% ylabel('Energy (total distances to the medoid in all clusters)')
% title("Multiparameter K-medoid cluster number vs performance of Kmedoid clustering with n=" + nwf + " waveform data")

