clear all; 
LoadCommonSettings_ericka;

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
NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor .* CS.SAR.data;
unNORMdata =  CS.SAR.data;

%% Apply SAR L1b_to_L2 code to find sigma0 

[DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi) ; 
sigma0 =  CS.sigma0.' ; 

%% Classify waveforms 
classification_method = 'Inger'; % choose from 'Inger', 'Rose', 'Peacock', ... 
class = Classify_Waveform(NORMdata, unNORMdata, sigma0, classification_method) ;  

%%Plot classification 
leadindex =  find(class == 2 ); 
iceindex = find(class == 1 ); 
ambindex = find(class == 0 ); 

figure
h= geoscatter( CS.GEO.LAT(iceindex),CS.GEO.LON(iceindex), 'filled', 'b') ;
hold on
geoscatter(CS.GEO.LAT(leadindex),CS.GEO.LON(leadindex),  'filled', 'r') ; 
%scatter( CS.GEO.LAT(ambindex),CS.GEO.LON(ambindex), 'filled', 'g')
hold off

