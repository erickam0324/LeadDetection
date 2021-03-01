% function LoadCommonSettings

% Load computer information to determine the path directories
[~, result] = system('hostname');
switch lower(deblank(result))
    
    case {'linux-ynpp.suse','cslobbelaptop.suse'}
        Drive = '/';
        PathHOME = fullfile(Drive, 'data', 'Projects', 'CryoSat2_Greenland', 'Software');
        PathDATA = fullfile(Drive, 'data', 'Projects', 'CryoSat2_Greenland', 'Data');
        PathDsof = fullfile(Drive, 'data', 'Projects', 'CryoSat2_Greenland', 'Software');       
        PathRSLT = fullfile(Drive, 'data', 'Projects', 'CryoSat2_Greenland', 'Results');
        Pth_FIG  = fullfile(Drive, 'data', 'Projects', 'CryoSat2_Greenland', 'Software');

    case {'tud278628'} % Ingers Desktop
        Drive = '/';
        PathHOME = fullfile(Drive, 'home', 'ibijdevaate', 'FAST4NL', 'Software');
        PathDATA = fullfile(Drive, 'data', 'FAST4NL', 'Data');
        PathDsof = fullfile(Drive, 'home', 'ibijdevaate', 'FAST4NL', 'Software');       
        PathRSLT = fullfile(Drive, 'data', 'FAST4NL', 'Results');
        Pth_FIG  = fullfile(Drive, 'home', 'ibijdevaate', 'FAST4NL', 'Software');  
        
    otherwise % Everything else
        % Check if it is windows or unix
        if isunix
            Drive = '/'; 
            PathHOME = fullfile(Drive, 'home', 'ibijdevaate', 'Projects', 'FAST4NL', 'Software');
            PathDATA = fullfile(Drive, 'home', 'ibijdevaate', 'Projects', 'FAST4NL', 'Data');
            PathDsof = fullfile(Drive, 'home', 'ibijdevaate', 'Projects', 'FAST4NL', 'Software');
            PathRSLT = fullfile(Drive, 'home', 'ibijdevaate', 'Projects', 'FAST4NL', 'Results');
            Pth_FIG  = fullfile(Drive, 'home', 'ibijdevaate', 'Projects', 'FAST4NL', 'Software');
        end 
end

%Specify path of work directory Matlab, here the m-files are stored
Pth_m     = fullfile(PathHOME, 'DART', filesep);

%Specify path where mat-files will be/are stored
Pth_mat   = fullfile(PathRSLT,'Matlab_Lib', filesep);

%Specify path where other software is stored
Pth_Soft  = fullfile(PathHOME, filesep);

%Specify path where GMT is located
Pth_GMT   = fullfile(PathRSLT,'GMT');

%Specify path where output of geopot07 is stored
Pth_geopot = fullfile(PathDsof, 'Geopot07', filesep);

if not(exist('Set','var'))
    Set = struct;
    Set.('MinLon')       = -17.5;
    Set.('MaxLon')       = 15.5;
    Set.('MinLat')       = 40.5;
    Set.('MaxLat')       = 66.5;
    Set.('Period')       = 365.2425;                %Period of seasonal cycle = 1 year
    Set.('PrintFigures') = true;
    Set.('PrintFormat')  = 'eps';
    Set.('PrintOpt')     = '';
    %Specify reference epoch (t0) that will be used to remove secular rates
    %in all kinds of data. These secular rates represent geophysical
    %phenomena like PGR, melting of the ice sheets, and sea level rise.
    Set.('REFdate')      = datenum([2013 07 01 00 00 00]);
    Set.('W0_EVRF2007')  = 62636857.89;             %Denker 2013 (Sciences of Geodesy - II: Innovations and Future Developments)
    Set.('W0_NAP')       = 62636857.25;             %Bursa et al. 2002 (World height system specified by geopotential at tide gauge stations)
    % Set.('W0_NAP')       = 62636857.02;             %Denker et al. 2004 (Status of the European Gravity and Geoid Project EGGP)
    Set.('W0_IERS2010')  = 62636856.0;              %Groten, 2004.
    Set.('RefEll')       = 'GRS80';

    Set.('ApplyRCR')      = true;                   %Apply Remove-Compute-Restore (RCR)
    Set.('RefFldLinrztn') = 'GOCO05s';              %Reference field used in linearization
    % Set.('RefFldLinrztn') = 'GO_CONS_GCF_2_DIR_R5'; %Reference field used in linearization
    % Set.('RefFldLinrztn') = 'DGM1S_19990101';       %Reference field used in linearization
    Set.('RCRfield')      = 'GOCO05s';              %Reference model used in RCR procedure
    % Set.('RCRfield')      = 'GO_CONS_GCF_2_DIR_R5'; %Reference model used in RCR procedure
    % Set.('RCRfield')      = 'DGM1S_19990101';       %Reference model used in RCR procedure
    Set.('Lmax_RCR')      = 280;                    %Set maximum degree of ref. model
    % Set.('Lmax_RCR')      = 300;                    %Set maximum degree of ref. model
    % Set.('Lmax_RCR')      = 250;                    %Set maximum degree of ref. model
    Set.('TS_RCR')        = 2;                      %Permanent tide system ref. model
end

[~,CallingFuncName,~] = star69;
if any(strcmp(CallingFuncName,{'PreProcessing_GRAVdata','EditShipboardDataCruisewise','BiasEstimation','Export2inputRBF','GetGRAVdataFrom','ApplyAlngTrSmoothingShipGrav','RTM_Settings'})) || strncmp(CallingFuncName,'GD_',3)
    %Names of airborne, terr., and/or shipboard gravity data providers
    GD_DatasetNames = {'RWS';'BGS_terr';'IFE';'BKG';'BGS_ship';'NKG';...
        'NGDC';'BGI_terr';'BGI_ship';'KSB';'GSNI';'DIAS';'NGI'};
    
    %Set #columns and column identifiers of global variable GravDATA
    colGD_NrCol      = 16;
    colGD_ProviderID = 1;  colGD_DataSetID  = 2;  colGD_TrackID    = 3;
    colGD_LAT        = 4;  colGD_LON        = 5;  colGD_h_P        = 6;
    colGD_h_Q        = 7;  colGD_h_Qast     = 8;  colGD_FAA        = 9;
    colGD_GDi        = 10; colGD_FAAref     = 11; colGD_GDiref     = 12; 
    colGD_ESTstd     = 13; colGD_TSAflag    = 14; colGD_RTM_FAA    = 15;
    colGD_RTM_GDi    = 16;
end

% Clear computer information
clear('result','Drive','CallingFuncName')
