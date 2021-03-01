%LoadCommonSettings              % In this script you should define your own folder structure, this will initialize variables such as PathDATA

%Computation settings
% defval('SAT','CS')                                                                %Satellite mission from which data are processed ('CS'/'S3A'/'S3B')
% defval('YY',[])                                                                   %Year for which you want to process data
% defval('MM',[])                                                                   %Month for which you want to process data
% defval('FNamesPreFix',[])                                                         %Set prefix of filenames that contain L1b data
% defval('DOM',[-90 90; -180 180])                                                  %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
% defval('Retracker','SAMOSA2')                                                     %Retracker to be used
% defval('SetRetr',{})

%Specify: 
%FName = '2017\04\CS_OFFL_SIR_SAR_1B_20170401T000352_20170401T000544_C001.DBL';                                              % filename of the dataset CS
%FName = '2017\04\S3A_SR_1_SRA____20170401T004834_20170401T013904_20170426T154541_3029_016_088______MAR_O_NT_002.SEN3';      % filename of the dataset CS
FName = 'S3A_SR_1_SRA____20201018T105409_20201018T110409_20201018T130551_0599_064_094______LN3_O_NR_004.SEN3';
% To read data
%[~,CS] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1',FName));       % This means your DBL files should be stored in 'PathDATA'/RadAlt/CryoSat/SIR_SAR_L1/.
[~,CS] = S3_L1b_read(FName)                                                             % within S3_L1b_read the full filepaths are defined as:
                                                                                                                                                                      
                                                                                        % PathL2  = fullfile(PathDATA,'RadAlt','Sentinel3A','SR_2_WAT');
                                                                                                                                                                          
% To process data                                                                                        
DATA = SAR_L1b_to_L2(SAT,FName,DOM,Retracker,SetRetr,1);
