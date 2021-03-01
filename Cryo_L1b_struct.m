function [HDR, CS]=Cryo_L1b_struct(full_filename,n_recs,DDAcf)

% Scope: Matlab Function to ingest CryoSat L1b data products to matlab workspace
% Data Level: 1a & 1b
% Supported Modes: FBR, LRM, SAR, SARin, FDM
%
% Input Argument: <full_filename> is the full pathname (path + file) where the CryoSat .DBL
% file is stored in your local drive
%
% Output Arguments: the structures <HDR>, containing the header of the read file, and the structure <CS>,
% containing the read data fields
%
% Author: Salvatore Dinardo
% Date: 14/04/2015
% Version: 1.4
% Compliance to CryoSat L1b Product Specification: 6.2
% Debugging: for any issues, please write to salvatore.dinardo@esa.int
% Track Change Log:
%
%    - version 1.1: fix for the fields:   CS.COR.TOTAL_gim and CS.COR.TOTAL_model
%    - version 1.3: add compatibility to Baseline C
%    - version 1.4:  nominal version for Baseline C

[pathname filename, ext]=fileparts(full_filename);
fid=fopen(fullfile(pathname,[filename ext ] ),'r','b');
s = dir(fullfile(pathname, [filename ext]  ));

%%%%%%%%%%%%%%%%%%%%%%%%  DATA HEADER READING %%%%%%%%%%%%%%%%%%%

MPH_size=1247;

i=0;
k=1;
while 1
    
    i=i+1;
    
    if ftell(fid)>MPH_size, break,   end
    tline = fgetl(fid);
    
    field=tline(1:strfind(tline,'=')-1);
    
    I=strfind(tline,'=')+1;
    
    if  i>2 && isfield(HDR,field)
        
        
        if strcmp(field,'DS_NAME')
            
            k=k+1;
            
        end
        
        field=[field num2str(k)];
        
    end
    
    if strcmp(tline(I),'"')
        
        value=tline(I+1:end-1);
        eval([ 'HDR.' field '='' ' value ''';']);
        
    else
        
        J=strfind(tline,'<')-1;
        
        if isempty(J)
            J=length(tline);
        end
        
        if  not(isempty(tline(I:J)))&& not(isnan(str2double(tline(I:J))))
            
            value=str2double(tline(I:J));
            eval([ 'HDR.' field '= ' num2str(value, '%10.5f') ';']);
            
        elseif not(isempty(tline(I:J)))&& (isnan(str2double(tline(I:J))))
            
            value=(tline(I:J));
            eval([ 'HDR.' field '= ''' value ''';']);
        end
        
    end
    
end

i=0;
k=1;

while 1
    
    i=i+1;
    
    if ftell(fid)>=MPH_size+HDR.SPH_SIZE
        break
    end
    
    tline = fgetl(fid);
    field=tline(1:strfind(tline,'=')-1);
    
    I=strfind(tline,'=')+1;
    
    if  i>2 && isfield(HDR,field)
        
        
        if strcmp(field,'DS_NAME')
            
            k=k+1;
            
        end
        
        field=[field num2str(k)];
        
    end
    
    if strcmp(tline(I),'"')
        
        value=tline(I+1:end-1);
        eval([ 'HDR.' field '='' ' value ''';']);
        
    else
        
        J=strfind(tline,'<')-1;
        if isempty(J)
            J=length(tline);
        end
        
        if  not(isempty(tline(I:J)))&& not(isnan(str2double(tline(I:J))))
            
            value=str2double(tline(I:J));
            eval([ 'HDR.' field '= ' num2str(value, '%10.5f') ';']);
            
        elseif not(isempty(tline(I:J)))&& (isnan(str2double(tline(I:J))))
            
            value=(tline(I:J));
            eval([ 'HDR.' field '= ''' value ''';']);
            
        end
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% OPERATIVE MODE IDENTIFICATION  %%%%%%%%%%%%%%%%%%

CS.GEO.OPERATION_MODE=strtrim(HDR.DS_NAME);
CS.GEO.Start_Time=(datenum(HDR.START_RECORD_TAI_TIME)-datenum('01-Jan-2000 00:00:00')).*24.*60.*60;
tmp=deblank(HDR.PRODUCT);
CS.GEO.Baseline_ID=tmp(end-3);

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INFORMATION %%%%%%%%%%%%%%%%%%%%%%

N_block=20;
SAR_pulses_burst=64;

switch  CS.GEO.Baseline_ID
    
    case  'C'
        
        N_samples=128;
%         N_samples_SAR=256;
        N_samples_SAR=N_samples*DDAcf.os_ZP;
        timestamp=4+4+4;
        time_group_FBR=timestamp+4+2+2+4*6+3.*4*3+4;  % Time and orbit group ( bytes)
        time_group_L1b=timestamp+4+2+2+4*6+3.*4*3+2+4*5;  % Time and orbit group ( bytes)
        GEO_extra_fields=2+4+4+4;
        
    otherwise
        
        N_samples=128;
        N_samples_SAR=128;
        timestamp=4+4+4;
        time_group_L1b=timestamp+4+2.*2+4.*6+3.*3.*4+4;  % Time and orbit group (84 bytes)
        time_group_FBR=timestamp+4+2.*2+4.*6+3.*3.*4+4;  % Time and orbit group (84 bytes)
        GEO_extra_fields=0;
        
end

measure_group=8+4.*18+4; % Measurement group (84 bytes)
geo_corr=4.*16;          % Geocorrection group (64 bytes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
        
        time_group= time_group_L1b;
        av_wfm_LRM=(4*6+8+2*N_samples+4*2+2*2); % Average Waveform group (300 bytes )
        wfm_LRM=2*N_samples+4*2+2*2; % Full Waveform group (268 bytes )
        record_size=geo_corr+av_wfm_LRM+(time_group+measure_group+wfm_LRM)*N_block;
        n_points=N_samples;
%         n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_FBR_SAR'
        
        GEO_extra_fields=0;
        time_group= time_group_FBR;
        wfm_SAR=2*N_samples*SAR_pulses_burst+2+2; % Waveform group (64 bytes )
        record_size=geo_corr+(time_group+measure_group+wfm_SAR)*N_block;
        n_points=2*N_samples*SAR_pulses_burst;
%         n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_FBR_SARIN'
        
        GEO_extra_fields=0;
        time_group= time_group_FBR;
        wfm_SIN=2.*2*(4.*N_samples)*SAR_pulses_burst+2+2;
        record_size=geo_corr+(time_group+measure_group+wfm_SIN)*N_block;
        n_points=2*(4.*N_samples)*SAR_pulses_burst;
%         n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_L1B_SAR'
        
        time_group= time_group_L1b;
        wfm_SAR=2*N_samples_SAR+4*2+2*2+50*2;
        av_wfm_SAR=(timestamp+4*3+8+2*N_samples+4*2+2*2);
        record_size=geo_corr+av_wfm_SAR+(time_group+measure_group+wfm_SAR)*N_block;
        n_points=N_samples_SAR;
%         n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    case 'SIR_L1B_SARIN'
        
        time_group= time_group_L1b;
        wfm_SIN=(N_samples_SAR*4)*2+4+4+2+2+50*2+(N_samples_SAR.*4).*2+(N_samples_SAR.*4).*4;
        av_wfm_SIN=12+4+4+4+8+N_samples.*4.*2+4+4+2+2;
        record_size=geo_corr+av_wfm_SIN+(time_group+measure_group+wfm_SIN)*N_block;
        n_points=N_samples_SAR.*4; % 512 bins insted of 128 bins
%         n_recs=(s.bytes-(MPH_size+HDR.SPH_SIZE))./record_size;
        
    otherwise
        
        disp('Mode not Supported or File Product not recognized');
        CS=[];HDR=[];return;
        
end

offset=geo_corr+(time_group+measure_group)*N_block;

if ~~(mod(n_recs,1))
    
    disp('File Product corrupt');
    HDR=[];CS=[];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%  DATA STRUCTURE INITIALIZATION %%%%%%%%%%%%%%%%%%%

CS.GEO.TAI.days=zeros(N_block,n_recs);
CS.GEO.TAI.secs=zeros(N_block,n_recs);
CS.GEO.TAI.microsecs=zeros(N_block,n_recs);
CS.GEO.USO=zeros(N_block,n_recs);
CS.GEO.MODE_ID=char(zeros(N_block,16,n_recs));
CS.GEO.SRC_CNT=zeros(N_block,n_recs);
CS.GEO.INS_CFG=char(zeros(N_block,32,n_recs));
CS.GEO.BURST_CNT=zeros(N_block,n_recs);
CS.GEO.LAT=zeros(N_block,n_recs);
CS.GEO.LON=zeros(N_block,n_recs);
CS.GEO.H=zeros(N_block,n_recs);
CS.GEO.H_rate=zeros(N_block,n_recs);
CS.GEO.V.Vx=zeros(N_block,n_recs);
CS.GEO.V.Vy=zeros(N_block,n_recs);
CS.GEO.V.Vz=zeros(N_block,n_recs);
CS.GEO.Beam.X=zeros(N_block,n_recs);
CS.GEO.Beam.Y=zeros(N_block,n_recs);
CS.GEO.Beam.Z=zeros(N_block,n_recs);
CS.GEO.BaseLine.X=zeros(N_block,n_recs);
CS.GEO.BaseLine.Y=zeros(N_block,n_recs);
CS.GEO.BaseLine.Z=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG=char(zeros(N_block,32,n_recs));


if  strcmp( CS.GEO.Baseline_ID, 'C')
    
    switch CS.GEO.OPERATION_MODE
        
        case {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
            
            CS.GEO.Star_Tracker_Usage=zeros(N_block,n_recs);
            CS.GEO.Antenna_Bench_Roll=zeros(N_block,n_recs);
            CS.GEO.Antenna_Bench_Pitch=zeros(N_block,n_recs);
            CS.GEO.Antenna_Bench_Yaw=zeros(N_block,n_recs);
            
    end
    
end

CS.GEO.MODE_ID_Tab.Ins_Oper_Mode=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Sarin_Degr_Case=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Reserved_1=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Cal4_Flag=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Plat_Att_Ctrl=zeros(N_block,n_recs);
CS.GEO.MODE_ID_Tab.Reserved_2=zeros(N_block,n_recs);


CS.GEO.INS_CFG_Tab.Rx_Chain_Use=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.SIRAL_ID=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_1=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Band_FLAG=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_2=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Tracking_Mode=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Ext_Cal=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Reserved_3=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Loop_Status=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Echo_Loss=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Real_Time_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Echo_Satur_Err=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Rx_Band_Atten=zeros(N_block,n_recs);
CS.GEO.INS_CFG_Tab.Cycle_Report_Err=zeros(N_block,n_recs);

switch  CS.GEO.Baseline_ID
    
    case  'C'
        
        CS.GEO.INS_CFG_Tab.Reserved_4=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Reserved_5=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Reserved_6=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.STR_ATTREF=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Reserved_7=zeros(N_block,n_recs);
        
    otherwise
        
        CS.GEO.INS_CFG_Tab.Star_Trk_1=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Star_Trk_2=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Star_Trk_3=zeros(N_block,n_recs);
        CS.GEO.INS_CFG_Tab.Reserved_4=zeros(N_block,n_recs);
        
end

CS.GEO.MCD_FLAG_Tab.Block_Degraded=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Blank_Block=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Datation_Degraded=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_Propag_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_File_Change=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Orbit_Discontinuity=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Echo_Saturation=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Other_Echo_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Rx1_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Rx2_Err_SARin=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Wind_Delay_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.AGC_Incon=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL1_Corr_Miss=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.DORIS_USO_Corr=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.Complex_CAL1_Corr_IPF=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.TRK_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.RX1_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.RX2_ECHO_Err=zeros(N_block,n_recs);
CS.GEO.MCD_FLAG_Tab.NPM_Incon=zeros(N_block,n_recs);

if strcmp(CS.GEO.Baseline_ID,'C')
    
    switch  CS.GEO.OPERATION_MODE
        
        case   {'SIR_L1B_LRM','SIR_L1B_FDM','SIR_L1B_SAR','SIR_L1B_SARIN'}
            
            CS.GEO.MCD_FLAG_Tab.CAL2_Corr_Missing=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_1=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_2=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_3=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Corr=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.CAL2_Corr_Miss=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.CAL2_Corr_IPF=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Power_Scaling_Err=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Att_Corr_Miss=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_4=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_5=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Mode=zeros(N_block,n_recs);
            
            
        case  {'SIR_FBR_SAR', 'SIR_FBR_SARIN'}
            
            CS.GEO.MCD_FLAG_Tab.Reserved_1=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Attitude_Corr_Missing=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.CAL1_Corr_Type=zeros(N_block,n_recs);
            CS.GEO.MCD_FLAG_Tab.Reserved_2=zeros(N_block,n_recs);
            
    end
    
else
    
    CS.GEO.MCD_FLAG_Tab.Reserved_1=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Reserved_2=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Reserved_3=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Reserved_4=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Corr=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.CAL2_Corr_Miss=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.CAL2_Corr_IPF=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Power_Scaling_Err=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Att_Corr_Miss=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Reserved_5=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Reserved_6=zeros(N_block,n_recs);
    CS.GEO.MCD_FLAG_Tab.Phase_Pertubation_Mode=zeros(N_block,n_recs);
    
end

CS.MEA.win_delay=zeros(N_block,n_recs);
CS.MEA.Ho=zeros(N_block,n_recs);
CS.MEA.Trk_H_rate=zeros(N_block,n_recs);
CS.MEA.LAI=zeros(N_block,n_recs);
CS.MEA.FAI=zeros(N_block,n_recs);
CS.MEA.AGC_1=zeros(N_block,n_recs);
CS.MEA.AGC_2=zeros(N_block,n_recs);
CS.MEA.Gain_Rx1=zeros(N_block,n_recs);
CS.MEA.Gain_Rx2=zeros(N_block,n_recs);
CS.MEA.Tx_Pwr=zeros(N_block,n_recs);
CS.MEA.dpl_range_corr=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_range_corr_rx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx_tx=zeros(N_block,n_recs);
CS.MEA.ins_gain_corr_rx=zeros(N_block,n_recs);
CS.MEA.int_phase_corr=zeros(N_block,n_recs);
CS.MEA.ext_phase_corr=zeros(N_block,n_recs);
CS.MEA.noise_power=zeros(N_block,n_recs);
CS.MEA.phase_slope_corr=zeros(N_block,n_recs);

CS.COR.dry_trop=zeros(1,n_recs);
CS.COR.wet_trop=zeros(1,n_recs);
CS.COR.inv_bar=zeros(1,n_recs);
CS.COR.dac=zeros(1,n_recs);
CS.COR.gim_ion=zeros(1,n_recs);
CS.COR.model_ion=zeros(1,n_recs);
CS.COR.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.ocean_loading_tide=zeros(1,n_recs);
CS.COR.solidearth_tide=zeros(1,n_recs);
CS.COR.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.surf_type=zeros(1,n_recs);
CS.COR.corr_status=char(zeros(32,n_recs));
CS.COR.corr_error=char(zeros(32,n_recs));


CS.COR.corr_status_tab.dry_trop=zeros(1,n_recs);
CS.COR.corr_status_tab.wet_trop=zeros(1,n_recs);
CS.COR.corr_status_tab.inv_bar=zeros(1,n_recs);
CS.COR.corr_status_tab.dac=zeros(1,n_recs);
CS.COR.corr_status_tab.gim_iono=zeros(1,n_recs);
CS.COR.corr_status_tab.model_iono=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_status_tab.surface_type=zeros(1,n_recs);
CS.COR.corr_status_tab.reserved=zeros(1,n_recs);

CS.COR.corr_error_tab.dry_trop=zeros(1,n_recs);
CS.COR.corr_error_tab.wet_trop=zeros(1,n_recs);
CS.COR.corr_error_tab.inv_bar=zeros(1,n_recs);
CS.COR.corr_error_tab.dac=zeros(1,n_recs);
CS.COR.corr_error_tab.gim_iono=zeros(1,n_recs);
CS.COR.corr_error_tab.model_iono=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_equilibrium_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_longperiod_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.ocean_loading_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.solidearth_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.geocentric_polar_tide=zeros(1,n_recs);
CS.COR.corr_error_tab.surface_type=zeros(1,n_recs);
CS.COR.corr_error_tab.reserved=zeros(1,n_recs);


switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
            
        CS.AVG.TAI.day=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(n_points,n_recs);
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.LRM.data=zeros(n_points.*N_block,n_recs);
        CS.LRM.echo_scaling=zeros(N_block,n_recs);
        CS.LRM.echo_scale_power=zeros(N_block,n_recs);
        CS.LRM.N_averaged_echoes=zeros(N_block,n_recs);
        CS.LRM.FLAG=char(zeros(16,N_block,n_recs));
        
    case {'SIR_L1B_SAR'}
        
        CS.AVG.TAI.days=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(N_samples,n_recs);
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.AVG.Mispointing_Err=zeros(1,n_recs);
        
        CS.SAR.data=zeros(n_points.*N_block,n_recs);
        CS.SAR.echo_scaling=zeros(N_block,n_recs);
        CS.SAR.echo_scale_power=zeros(N_block,n_recs);
        CS.SAR.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SAR.FLAG=char(zeros(16,N_block,n_recs));
        CS.SAR.FLAG_tab.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SAR.FLAG_tab.Reserved=zeros(N_block,n_recs);
        CS.SAR.beam_param=zeros(50,N_block,n_recs);
            
    case {'SIR_FBR_SAR'}
        
        CS.FBR.data=zeros(n_points.*N_block,n_recs,'int8');
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case {'SIR_FBR_SARIN'}
        
        CS.FBR.data_ch1=zeros(n_points.*N_block,n_recs,'int8');
        CS.FBR.data_ch2=zeros(n_points.*N_block,n_recs,'int8');
        CS.FBR.N_pulses=zeros(N_block,n_recs);
        CS.FBR.FLAG=char(zeros(16,N_block,n_recs));
        
    case 'SIR_L1B_SARIN'
        
        CS.AVG.TAI.days=zeros(1,n_recs);
        CS.AVG.TAI.secs=zeros(1,n_recs);
        CS.AVG.TAI.microsecs=zeros(1,n_recs);
        CS.AVG.lat=zeros(1,n_recs);
        CS.AVG.lon=zeros(1,n_recs);
        CS.AVG.H=zeros(1,n_recs);
        CS.AVG.win_delay=zeros(1,n_recs);
        CS.AVG.data=zeros(N_samples.*4,n_recs);
        CS.AVG.echo_scaling=zeros(1,n_recs);
        CS.AVG.echo_scale_power=zeros(1,n_recs);
        CS.AVG.N_averaged_echoes=zeros(1,n_recs);
        CS.AVG.OneHz_Echo_Err=zeros(1,n_recs);
        CS.AVG.Mispointing_Err=zeros(1,n_recs);
        
        CS.SIN.data=zeros(n_points.*N_block,n_recs);
        CS.SIN.echo_scaling=zeros(N_block,n_recs);
        CS.SIN.echo_scale_power=zeros(N_block,n_recs);
        CS.SIN.N_averaged_echoes=zeros(N_block,n_recs);
        CS.SIN.FLAG=char(zeros(16,N_block,n_recs));
        CS.SIN.FLAG_tab.Approximate_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Exact_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Doppler_Weighting_Computed=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Doppler_Weighting_Applied_Before_Stack=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Multilook_Incomplete=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Beam_Angle_Steering_Err=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.AntiAliased_Power_Echo=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Auto_Beam_Steering=zeros(N_block,n_recs);
        CS.SIN.FLAG_tab.Reserved=zeros(N_block,n_recs);
        CS.SIN.beam_param=zeros(50,N_block,n_recs);
        CS.SIN.coherence=zeros(n_points.*N_block,n_recs);
        CS.SIN.phase_difference=zeros(n_points.*N_block,n_recs);
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CS.GEO.V.V=sqrt(CS.GEO.V.Vx.^2+CS.GEO.V.Vy.^2+CS.GEO.V.Vz.^2);
CS.GEO.Serial_Sec_Num=CS.GEO.TAI.days.*24.*60.*60+CS.GEO.TAI.secs+CS.GEO.TAI.microsecs./1e6;

CS.COR.TOTAL_gim=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.dac+CS.COR.gim_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
    CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.COR.TOTAL_model=CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.dac+CS.COR.model_ion+CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+...
    CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide;

CS.GEO.Elapsed_Time=zeros(size(CS.GEO.Serial_Sec_Num));
CS.GEO.Elapsed_Time(CS.GEO.Serial_Sec_Num~=0)=CS.GEO.Serial_Sec_Num(CS.GEO.Serial_Sec_Num~=0)-CS.GEO.Start_Time;

CS.GEO.MODE_ID=CS.GEO.MODE_ID_Tab;
CS.GEO=rmfield(CS.GEO,'MODE_ID_Tab');

CS.GEO.INS_CFG=CS.GEO.INS_CFG_Tab;
CS.GEO=rmfield(CS.GEO,'INS_CFG_Tab');

CS.GEO.MCD_FLAG=CS.GEO.MCD_FLAG_Tab;
CS.GEO=rmfield(CS.GEO,'MCD_FLAG_Tab');

CS.COR.corr_status=CS.COR.corr_status_tab;
CS.COR=rmfield(CS.COR,'corr_status_tab');

CS.COR.corr_error=CS.COR.corr_error_tab;
CS.COR=rmfield(CS.COR,'corr_error_tab');

switch CS.GEO.OPERATION_MODE
    
    case {'SIR_L1B_LRM','SIR_L1B_FDM'}
        
        CS.LRM.data=reshape(CS.LRM.data,N_samples,N_block,n_recs);
        
    case {'SIR_L1B_SAR'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
        CS.SAR.data=reshape(CS.SAR.data,n_points,N_block,n_recs);
        
        CS.SAR.FLAG=CS.SAR.FLAG_tab;
        CS.SAR=rmfield(CS.SAR,'FLAG_tab');
        
    case {'SIR_L1B_SARIN'}
        
        CS.AVG.Serial_Sec_Num=CS.AVG.TAI.days.*24.*60.*60+CS.AVG.TAI.secs+CS.AVG.TAI.microsecs./1e6;
        CS.AVG.Elapsed_Time=CS.AVG.Serial_Sec_Num-CS.GEO.Start_Time;
        CS.SIN.data=reshape(CS.SIN.data,n_points,N_block,n_recs);
        CS.SIN.coherence=reshape(CS.SIN.coherence,n_points,N_block,n_recs);
        CS.SIN.phase_difference=reshape(CS.SIN.phase_difference,n_points,N_block,n_recs);
        CS.SIN.data=reshape(CS.SIN.data,n_points,N_block,n_recs);
        
        CS.SIN.FLAG=CS.SIN.FLAG_tab;
        CS.SIN=rmfield(CS.SIN,'FLAG_tab');
        
    case  {'SIR_FBR_SAR'}
        
        CS.FBR.data=reshape((complex(CS.FBR.data(1:2:end,:),CS.FBR.data(2:2:end,:))),N_samples,SAR_pulses_burst,N_block,[]);
        
    case  {'SIR_FBR_SARIN'}
        
        CS.FBR.data_ch1=reshape((complex(CS.FBR.data_ch1(1:2:end,:),CS.FBR.data_ch1(2:2:end,:))),4.*N_samples,SAR_pulses_burst,[]);
        CS.FBR.data_ch2=reshape((complex(CS.FBR.data_ch2(1:2:end,:),CS.FBR.data_ch2(2:2:end,:))),4.*N_samples,SAR_pulses_burst,N_block,[]);
        
end

fclose(fid);
