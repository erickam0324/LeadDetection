function [DATA,CS] = SAR_L1b_to_L2_ericka(SAT,FName,DOM,Retracker,SetRetr,DEM,IDXpoi)

%SAR_L1B_TO_L2 processes level 2 data from CryoSat & Sentinel-3A/B level 1b
%SAR data.

%% Settings
%General settings
%LoadCommonSettings

%Computation settings
defval('SAT','CS')                                                                %Satellite mission from which data are processed ('CS'/'S3A'/'S3B')
defval('FName','2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001.DBL') %*.DBL/*.nc file that contains level 1b data
defval('DOM',[60 90; -180 180])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','SAMOSA2')                                                   %Retracker to be used
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %s;',SetRetr{i,1},SetRetr{i,2})); end
defval('SAMOSAimpl','SD')                                                         %SAMOSA implementation to be used ('SD' = one described by Salvatore Dinardo et al. 2018; 'DPM_v2_5_2' = Gommenginger et al. 2017)
defval('SolverAlgo','TRR')                                                        %Solver Algorithm to be applied in case Retracker is an analytical retracker ('LM' = Levenberg-Marquardt; 'TRR' = trust-region-reflective)
defval('max_waveform_energy',150)                                                 %Maximum allowed energy of NORMALIZED waveform (i.e., area below the curve)
defval('ApplyWfClass',true)                                                       %Apply surface type classification based on waveform parameters [true/false].
defval('MAfromSTR',true)                                                          %Obtain mispointing angle from star tracker data
defval('MINbin',10);                                                              %Minimum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('MAXbin',125);                                                             %Maximum bin index of interval in which retracking point is accepted for NON-zero padded waveform
defval('Wf_Norm_Aw',1)                                                            %Width (in gates) of the sliding window for waveform normalization
defval('Wf_TN_First',3);                                                          %First gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('Wf_TN_Last',8);                                                           %Last gate of the NON-zero padded waveform to estimate the amplitude of the thermal noise floor
defval('UseDEM2CompInitRP',false)                                                 %Use DEM to compute initial retracking point [true/false].
defval('DEM',[])                                                                  %DEM used to compute initial retracking point
if UseDEM2CompInitRP && isempty(DEM), error('Provide DEM as input!'), end

%Create default optimization options for SolverAlgo with the named
%parameters altered with the specified values
if ~any(strcmp(Retracker,{'OCOG','Threshold'}))
    switch SolverAlgo
        case 'LM'
            %Levenberg-marquardt
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt','StepTolerance',1E-2,'FunctionTolerance',1E-2,'OptimalityTolerance',1E-2);
        case 'TRR'
            %Trust-region-reflective
            options = optimoptions('lsqcurvefit','Display','off','Algorithm','trust-region-reflective','StepTolerance',1E-2,'FunctionTolerance',1E-2,'OptimalityTolerance',1E-2);
        otherwise
            error('Solver not reckognized!')
    end
end

%Physical constants & instrument characteristics
DDAcf   = DDA_ConfigFile(SAT,'SAR');

%Remaining settings
defval('MakePlots',exist('IDXpoi','var')) %Make plots
if MakePlots; if ~isscalar(IDXpoi); MakePlots = false; end; end
defval('FntSz',14)                        %Set fontsize to be used in figures

%% Read & Crop Level 1b SAR data
switch SAT
    case {'CS','CryoSat'}
        %Read data
        if isstruct(FName)
            CS     = FName;
        else
            [~,CS] = Cryo_L1b_read(fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1',FName));
        end

        %Wrap angle in degrees to [-180 180]
        CS.GEO.LON = wrapTo180(CS.GEO.LON);
        %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
        if max(abs(diff(CS.GEO.LON(CS.SAR.N_averaged_echoes(:) > 0)))) > 350
            CS.GEO.LON = wrapTo360(CS.GEO.LON);
        end
        
        %Remove field AVG to save RAM
        CS = rmfield(CS,'AVG');
        
        %If it does not exists, set field with beam angles
        if ~isfield(CS.SAR,'BeamAngle')
            CS.SAR.('BeamAngle') = cell(numel(CS.GEO.LAT),1);
        end
        
        %Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
        CS.MEA.win_delay = CS.MEA.win_delay.*(CS.GEO.USO+1);
        
        %Compute reference range
        CS.MEA.('ref_range') = 0.5*DDAcf.c*CS.MEA.win_delay(:);
        
        %Transform acquisition time to datenum format
        CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400 + CS.GEO.TAI.microsecs(:)./1e6./86400 - DDAcf.Biases.Datation/86400;
        
        %Get Roll, Yaw, and Pitch angles from star tracker data
        if MAfromSTR
            try
                tSTR          = datetime(CS.GEO.Start_Time./(24.*60.*60) + datenum('01-Jan-2000 00:00:00') - 1,'ConvertFrom','datenum');
                mpSTR         = READ_STR_Mispointing_Angles(tSTR);
                DUM           = interp1(mpSTR.Time,[mpSTR.Pitch mpSTR.Roll mpSTR.Yaw],CS.TIME+DDAcf.Biases.Datation/86400,'spline');
                IDXnan        = isnan(interp1(mpSTR.Time,mpSTR.Roll,CS.TIME+DDAcf.Biases.Datation/86400));
                DUM(IDXnan,:) = NaN;
                CS.GEO.Antenna_Bench_Roll(~isnan(DUM(:,2)))  = DUM(~isnan(DUM(:,2)),2);
                %Note that yaw angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Yaw(~isnan(DUM(:,3)))   = -DUM(~isnan(DUM(:,3)),3);
                %Note that pitch angle is defined in different sign convention!
                CS.GEO.Antenna_Bench_Pitch(~isnan(DUM(:,1))) = -DUM(~isnan(DUM(:,1)),1);
            catch exception
                fprintf('%s\n',exception.message)
            end
        end
        
        % %Salvatore Dinardo applies the pitch and roll biases proposed by Remko
        % %Scharroo. In doing so, he undoes the biases applied in producing baseline
        % %C (Main evolutions and expected quality improvements in Baseline C Level
        % %1b products - ARESYS / ESA, V1.3, https://wiki.services.eoportal.org/tiki-
        % %index.php?page=CryoSat+Technical+Notes)
        % CS.GEO.Antenna_Bench_Roll  = CS.GEO.Antenna_Bench_Roll  + 0.1062 - 0.086;
        % CS.GEO.Antenna_Bench_Pitch = CS.GEO.Antenna_Bench_Pitch + 0.0550 - 0.096;
        
        %Convert Roll, Yaw & Pitch angles to radians
        CS.GEO.Antenna_Bench_Roll  = deg2rad(CS.GEO.Antenna_Bench_Roll);
        CS.GEO.Antenna_Bench_Yaw   = deg2rad(CS.GEO.Antenna_Bench_Yaw);
        CS.GEO.Antenna_Bench_Pitch = deg2rad(CS.GEO.Antenna_Bench_Pitch);
        
        %The power echo	sample values are all scaled to	fit	between	0 and 65535.
        %The scaling factors can change	for	each waveform. To convert these back to
        %values in Watts the following equation should be used (CryoSat Product
        %Handbook, Eqn 4.2�?1):
        %Power in Watts	= scaled value * (scale	factor * 10^�?9) * 2^scale power
        CS.SAR.data = bsxfun(@times,CS.SAR.data,reshape((CS.SAR.echo_scaling.*1E-9) .* 2.^(CS.SAR.echo_scale_power),1,size(CS.SAR.data,2),size(CS.SAR.data,3)));
        
        %Identify flagged data
        IDXfd = CS.GEO.MCD_FLAG.Block_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Blank_Block(:) == 1 | ...
            CS.GEO.MCD_FLAG.Datation_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Orbit_Propag_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.Echo_Saturation(:) == 1 | CS.GEO.MCD_FLAG.Other_Echo_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.Rx1_Err_SARin(:) == 1 | CS.GEO.MCD_FLAG.Rx2_Err_SARin(:) == 1 | ...
            CS.GEO.MCD_FLAG.Wind_Delay_Incon(:) == 1 | CS.GEO.MCD_FLAG.AGC_Incon(:) == 1 | ...
            CS.GEO.MCD_FLAG.TRK_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.RX1_ECHO_Err(:) == 1 | ...
            CS.GEO.MCD_FLAG.RX2_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.NPM_Incon(:) == 1 | ...
            CS.GEO.MCD_FLAG.Power_Scaling_Err(:) == 1;
        
        %Compute sum of instrumental, propagation and geophysical corrections
        %Corrections to be applied in case (i) surface == open oceans or
        %semi�?enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
        %still lacking). Results in SSH as defined in Eq. 1 (Dinardo et al. 2018)
        CorrST01  = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
            CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion+...
            CS.COR.dac;
        %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
        %Dinardo et al. 2018) is needed. Note that the SSB correction still has to
        %be applied.
        CorrST01i = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+CS.COR.gim_ion+...
            CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide;
        %Corrections to be applied in case (i) surface == continental ice OR (ii)
        %surface == land
        CorrST23  = CS.COR.ocean_loading_tide+CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+...
            CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.gim_ion;
        warning('In case of sea ice, different set of corrections should be applied!')

        %Interpolate corrections to 20 Hz
        CS.('surf_type') = reshape(repmat(CS.COR.surf_type,20,1),numel(CS.COR.surf_type)*20,1);
        CorrST01         = reshape(repmat(CorrST01,20,1),numel(CorrST01)*20,1);
        CorrST01i        = reshape(repmat(CorrST01i,20,1),numel(CorrST01i)*20,1);
        CorrST23         = reshape(repmat(CorrST23,20,1),numel(CorrST23)*20,1);
%         CS.COR.SAL       = reshape(repmat(CS.COR.ocean_loading_tide,20,1),numel(CS.COR.ocean_loading_tide)*20,1);
        FlgC01           = isnan(CorrST01); FlgC01i = isnan(CorrST01i); FlgC23 = isnan(CorrST23);
        cEqTide          = reshape(repmat(CS.COR.ocean_equilibrium_tide,20,1),numel(CS.COR.ocean_equilibrium_tide)*20,1);
        cLPTide          = reshape(repmat(CS.COR.ocean_longperiod_tide,20,1),numel(CS.COR.ocean_longperiod_tide)*20,1);
        cOLTide          = reshape(repmat(CS.COR.ocean_loading_tide,20,1),numel(CS.COR.ocean_loading_tide)*20,1);
        cSETide          = reshape(repmat(CS.COR.solidearth_tide,20,1),numel(CS.COR.solidearth_tide)*20,1);
        cGPTide          = reshape(repmat(CS.COR.geocentric_polar_tide,20,1),numel(CS.COR.geocentric_polar_tide)*20,1);
        cMeteo           = reshape(repmat(CS.COR.dac,20,1),numel(CS.COR.dac)*20,1);      
        
	CS.COR.LAT1Hz = CS.GEO.LAT(1,:);
	CS.COR.LON1Hz = CS.GEO.LON(1,:);
  
    case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
        %Read data
        if isstruct(FName)
            CS     = FName;
        else
            [~,CS] = S3_L1b_read_ericka(FName);
            if isempty(CS), DATA.('LAT') = []; return; end
        end

        %Wrap angle in degrees to [-180 180]
        CS.GEO.LON = wrapTo180(CS.GEO.LON);
        %If segment crosses -180 meridian, wrap angle in degrees to [0 360]
        if max(abs(diff(CS.GEO.LON(CS.SAR.N_averaged_echoes(:) > 0)))) > 350
            CS.GEO.LON = wrapTo360(CS.GEO.LON);
        end
        
        %Convert matrix with beam angles to cell array
        CS.SAR.BeamAngle = mat2cell(CS.SAR.BeamAngle, size(CS.SAR.BeamAngle,1), ones(1,size(CS.SAR.BeamAngle,2)));
        
        %Compute tangential velocities
        CS.GEO.V.('V') = sqrt(CS.GEO.V.Vx.^2 + CS.GEO.V.Vy.^2 + CS.GEO.V.Vz.^2);

        %Transform acquisition time to datenum format
        CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400;

        %Set Roll, Yaw & Pitch angles to radians
        CS.GEO.Antenna_Bench_Roll  = zeros(size(CS.GEO.LAT));
        CS.GEO.Antenna_Bench_Yaw   = zeros(size(CS.GEO.LAT));
        CS.GEO.Antenna_Bench_Pitch = zeros(size(CS.GEO.LAT));

        %The echo is corrected for Doppler range effect, phase/power burst
        %calibration and GPRW effect . the echo is scaled using the
        %corrected AGC (agc_ku_l1b_echo_sar_ku). Note that ''it is not
        %possible to convert power waveform in watts because it requires
        %some instrumental information from the industry who designed the
        %altimeter and this information is not available to the users.
        %Nevertheless, the S3A waveforms can be compared to each other, by
        %applying the agc values provided in the SRAL L1 products.''

        %Identify flagged data
        IDXfd = CS.GEO.FLAGS.time_corr_val_l1b_echo_sar_ku(:) == 2 | ...   %time correlation info invalid % CS.GEO.FLAGS.man_pres_l1b_echo_sar_ku(:) ~= 0 | ...            %ongoing manoeuvre
            CS.GEO.FLAGS.man_thrust_l1b_echo_sar_ku(:) == 1 | ...          %ongoing thrust % CS.GEO.FLAGS.man_plane_l1b_echo_sar_ku(:) == 1 | ...           %out of plane
            CS.GEO.FLAGS.gnss_status_l1b_echo_sar_ku(:) == 1 | ...         %navigation message GNSS receiver not valid/available
            CS.GEO.FLAGS.loss_track_l1b_echo_sar_ku(:) == 1;               %loss of track

	% Replace missing GIM iono corrections by altimeter derived
	CS.COR.iono_cor_gim_01_ku(isnan(CS.COR.iono_cor_gim_01_ku)) = CS.COR.iono_cor_alt_01_ku(isnan(CS.COR.iono_cor_gim_01_ku));

        %Compute sum of instrumental, propagation and geophysical corrections
        %Corrections to be applied in case (i) surface == open oceans or
        %semi�?enclosed seas OR (ii) surface == enclosed seas or lakes.
        %Results in SSH as defined in Eq. 1 (Dinardo et al. 2018). Note
        %that the SSB correction still has to be applied as the one that is
        %in the L2 file has not been tuned for S3A/S3B and contains Jason-2
        %SSB solution (Sentinel-3 Product Notice – STM L2 Marine, 6
        %February 2019, Notice #5)
        CorrST01           = CS.COR.ocean_equilibrium_tide+CS.COR.ocean_longperiod_tide+CS.COR.ocean_loading_tide+...
            CS.COR.solidearth_tide+CS.COR.geocentric_polar_tide+CS.COR.dry_trop+CS.COR.wet_trop+...
            CS.COR.inv_bar+CS.COR.hf_fluct_cor+CS.COR.iono_cor_gim_01_ku+...
            CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;
        IDXnn              = knnsearch(CS.COR.TIME1Hz,CS.TIME);
        FlgC01             = abs(CS.TIME-CS.COR.TIME1Hz(IDXnn)) > 1/86400 | isnan(CorrST01(IDXnn)) | isnan(CS.COR.LAT1Hz(IDXnn));
        
        %Corrections to be applied in case the instantaneous SSHs (SSHi, see Eq. 4
        %Dinardo et al. 2018) are needed. Note that the SSB correction
        %still has to be applied as the one that is in the L2 file has not
        %been tuned for S3A/S3B and contains Jason-2 SSB solution
        %(Sentinel-3 Product Notice – STM L2 Marine, 6 February 2019,
        %Notice #5)
        CorrST01i          = CS.COR.dry_trop+CS.COR.wet_trop+CS.COR.ocean_loading_tide+...
            CS.COR.solidearth_tide+0.468*CS.COR.geocentric_polar_tide+CS.COR.iono_cor_gim_01_ku+...
            CS.COR.cog_cor_01+CS.COR.mod_instr_cor_range_01_ku;
        FlgC01i            = abs(CS.TIME-CS.COR.TIME1Hz(IDXnn)) > 1/86400 | isnan(CorrST01i(IDXnn)) | isnan(CS.COR.LAT1Hz(IDXnn));
        
        
        if nnz(~isnan(CorrST01)) >= 2 && nnz(~isnan(CorrST01i)) >= 2   
            try 
                CorrST01       = interp1(CS.COR.TIME1Hz(~isnan(CorrST01)),CorrST01(~isnan(CorrST01)),CS.TIME,'pchip',NaN);
                CorrST01i      = interp1(CS.COR.TIME1Hz(~isnan(CorrST01i)),CorrST01i(~isnan(CorrST01i)),CS.TIME,'pchip',NaN);           
                
                cEqTide        = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.ocean_equilibrium_tide)),CS.COR.ocean_equilibrium_tide(~isnan(CS.COR.ocean_equilibrium_tide)),CS.TIME,'pchip',NaN);
                cLPTide        = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.ocean_longperiod_tide)),CS.COR.ocean_longperiod_tide(~isnan(CS.COR.ocean_longperiod_tide)),CS.TIME,'pchip',NaN);
                cOLTide        = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.ocean_loading_tide)),CS.COR.ocean_loading_tide(~isnan(CS.COR.ocean_loading_tide)),CS.TIME,'pchip',NaN);
                cSETide        = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.solidearth_tide)),CS.COR.solidearth_tide(~isnan(CS.COR.solidearth_tide)),CS.TIME,'pchip',NaN);
                cGPTide        = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.geocentric_polar_tide)),CS.COR.geocentric_polar_tide(~isnan(CS.COR.geocentric_polar_tide)),CS.TIME,'pchip',NaN);
                cMeteo         = interp1(CS.COR.TIME1Hz(~isnan(CS.COR.inv_bar)),CS.COR.inv_bar(~isnan(CS.COR.inv_bar)),CS.TIME,'pchip',NaN);
            catch
                [~,u,~] = unique(CS.COR.TIME1Hz);
                if length(u) < 3
                    return
                end
                CorrST01       = interp1(CS.COR.TIME1Hz(u),CorrST01(u),CS.TIME,'pchip',NaN);
                CorrST01i      = interp1(CS.COR.TIME1Hz(u),CorrST01i(u),CS.TIME,'pchip',NaN);           
                               
                
                cEqTide        = interp1(CS.COR.TIME1Hz(u),CS.COR.ocean_equilibrium_tide(u),CS.TIME,'pchip',NaN);
                cLPTide        = interp1(CS.COR.TIME1Hz(u),CS.COR.ocean_longperiod_tide(u),CS.TIME,'pchip',NaN);
                cOLTide        = interp1(CS.COR.TIME1Hz(u),CS.COR.ocean_loading_tide(u),CS.TIME,'pchip',NaN);
                cSETide        = interp1(CS.COR.TIME1Hz(u),CS.COR.solidearth_tide(u),CS.TIME,'pchip',NaN);
                cGPTide        = interp1(CS.COR.TIME1Hz(u),CS.COR.geocentric_polar_tide(u),CS.TIME,'pchip',NaN);
                cMeteo         = interp1(CS.COR.TIME1Hz(u),CS.COR.inv_bar(u),CS.TIME,'pchip',NaN);

            end
            CS.COR.SAL 	   = zeros(size(CS.TIME));
        else
            CorrST01       = zeros(size(CS.TIME));
            CorrST01i      = zeros(size(CS.TIME));
            cEqTide        = zeros(size(CS.TIME));
            cLPTide        = zeros(size(CS.TIME));
            cOLTide        = zeros(size(CS.TIME));
            cSETide        = zeros(size(CS.TIME));
            cGPTide        = zeros(size(CS.TIME));
            cMeteo         = zeros(size(CS.TIME));
       	    CS.COR.SAL 	   = zeros(size(CS.TIME));
        end
%         CorrST01(FlgC01)   = 0;

%         if nnz(~isnan(CorrST01i)) >= 2
%             try
%                 CorrST01i      = interp1(CS.COR.TIME1Hz(~isnan(CorrST01i)),CorrST01i(~isnan(CorrST01i)),CS.TIME,'pchip',NaN);
%             catch
%                 [~,u,~] = unique(CS.COR.TIME1Hz);
%                 CS.COR.TIME1Hz = CS.COR.TIME1Hz(u);
%                 CorrST01i = CorrST01i(u);
%                 CorrST01i      = interp1(CS.COR.TIME1Hz(~isnan(CorrST01i)),CorrST01i(~isnan(CorrST01i)),CS.TIME,'pchip',NaN);
%             end
%         else
%             CorrST01i      = zeros(size(CS.TIME));
%         end
%         CorrST01i(FlgC01i) = 0;
        %Corrections to be applied in case (i) surface == continental ice OR (ii)
        %surface == land
        CorrST23           = zeros(size(CS.GEO.LAT));
        FlgC23             = true(size(CS.GEO.LAT));
        warning('In case of sea ice, different set of corrections should be applied!')

    otherwise
        error('SAT: %s not recognized',SAT)
end

%Determine Zero-Padding Oversampling Factor and Waveform sampling interval.
%If os_ZP ~= 1, adjust reference bin.
NrBins = size(CS.SAR.data,1);        %Nr of bins/samples in any waveform
os_ZP  = NrBins/DDAcf.Np;            %Zero-Padding Oversampling Factor
RefBin = (DDAcf.RefBin-1)*os_ZP + 1;
dt     = 1/(os_ZP*DDAcf.B);          %Waveform sampling interval [s]

% %Compute mispointing angles
% CS.GEO.('MPA') = acos(cos(CS.GEO.Antenna_Bench_Pitch)./sqrt(sum(cat(3,sin(CS.GEO.Antenna_Bench_Roll),sin(CS.GEO.Antenna_Bench_Yaw),cos(CS.GEO.Antenna_Bench_Pitch)).^2,3)));

%Compute the local radii of curvature of the Earth's surface (Maulik Jain,
%Improved sea level determination in the Arctic regions through development
%of tolerant altimetry retracking, Eq. 5.2, pp. 47)
CS.GEO.('Re')  = sqrt(DDAcf.RefEll.SemimajorAxis^2*cosd(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2 + DDAcf.RefEll.SemiminorAxis^2*sind(geocentricLatitude(CS.GEO.LAT,DDAcf.RefEll.Flattening)).^2);

%Compute slope of orbit
if CS.GEO.LAT(find(~IDXfd,1,'first')) < CS.GEO.LAT(find(~IDXfd,1,'last'))
    track_sign = -1; %ascending track
else
    track_sign = 1;  %descending track
end
CS.GEO.('orbit_slope') = track_sign*((DDAcf.RefEll.SemimajorAxis^2 - DDAcf.RefEll.SemiminorAxis^2)./(2*CS.GEO.Re.^2)).*sin(2*deg2rad(CS.GEO.LAT)) - (-CS.GEO.H_rate./CS.GEO.V.V);

%Normalize each waveform by max value of waveform
NORMfactor = 1./max(movmean(CS.SAR.data,Wf_Norm_Aw,1),[],1);
NORMdata   = NORMfactor .* CS.SAR.data;

%Estimate the normalized thermal noise
CS.MEA.('TN') = mean(NORMdata(Wf_TN_First*os_ZP:Wf_TN_Last*os_ZP,:,:));

%Crop data to area of interest
if ischar(DOM)
    switch DOM
        case {'JakobsHavn','jakobshavn'}
            %Load polygon that outlines the JakobsHavn glacier
            load(fullfile(PathDATA,'Geography','jakobshavn_basin.mat'))
            IDX = inpolygon(CS.GEO.LON(:),CS.GEO.LAT(:),polyLon,polyLat);
        otherwise
            error('DOMain of interest not reckognized')
    end
else
    IDX = ingeoquad(CS.GEO.LAT(:),CS.GEO.LON(:),DOM(1,:),DOM(2,:));
end

%% Data editing
%Exclude flagged data
IDX(IDXfd) = false;

%Identify waveforms for which power == 0 for all entries
IDX(squeeze(all(CS.SAR.data == 0,1))) = false;

%Energy NORMALIZED waveform < max_waveform_energy && power at first bins of
%NORMALIZED waveform should be at noise level
IDX(squeeze(trapz(1:NrBins,NORMdata,1) >= max_waveform_energy) | squeeze(any(NORMdata(1:MINbin*os_ZP,:,:) > .1))) = false;

%Return if no data remain
if ~any(IDX)
    DATA = struct('TIME',[],'LAT',[],'LON',[],'HEI',[],'SurfT',[]);
    return
end

%Select points of interest (by default entire track is retracked)
defval('IDXpoi',find(IDX)')

% %% Classify waveforms
% if ApplyWfClass
%     %Preliminaries
     [CS.('n_ret'),CS.('Pu'),CS.('WFc')] = deal(nan(numel(CS.GEO.LAT),1));
    
    %Apply threshold retracker to obtain retracking point and Pu
    [CS.n_ret(IDXpoi),CS.Pu(IDXpoi)]  = RETRACKER.Threshold_mat(NORMdata(:,IDXpoi));
    
    %Compute range and backscatter coefficient (sigma0 / sigma-naught)
    CS.('range')   = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
    switch SAT
        case {'CS','CryoSat'}
            CS.Pu          = CS.Pu ./ NORMfactor(:);
            CS.('sigma0')  = Compute_Sigma0(CS,DDAcf);
            CS.Pu          = CS.Pu .* NORMfactor(:);
        case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
            %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
            %Definition, Accuracy and Specification, Section 2.15.3.3
             CS.('sigma0')  = CS.SAR.scale_factor_ku_l1b_echo_sar_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
        otherwise
            error('SAT: %s not recognized',SAT)
    end
    
%     %Classify waveforms
%     CS.WFc(IDXpoi) = Classify_Waveforms(SAT,NORMdata(:,IDXpoi),CS.sigma0(IDXpoi)',CS.surf_type(IDXpoi)');
% else
%     CS.WFc(IDXpoi) = CS.surf_type(IDXpoi);
% end
%     
% %% Preliminaries
% %Declare arrays
% [CS.('n_ret'),CS.('Pu'),CS.('SWH')] = deal(nan(numel(CS.GEO.LAT),1));
% [CS.('nu'),CS.('ExitF'),CS.('MF')]  = deal(nan(numel(CS.GEO.LAT),1));
% CS.('PCorr')                        = nan(numel(CS.GEO.LAT),1);
% BinIDs                              = (1:NrBins)';
% CS.('WDrecon')                      = nan(size(CS.SAR.data));
% 
% %Define/load look-up-tables
% if strncmp(Retracker,'SAMOSA',6)
%     if strcmp(SAMOSAimpl,'SD')
%         %Generate look-up-tables for fast evaluation of modified Bessel
%         %functions of the first kind
%         LUT_x    = logspace(-16,4,10000)';
%         LUT_B14  = griddedInterpolant(LUT_x,besseli(1/4,LUT_x,1),'spline');
%         LUT_Bm14 = griddedInterpolant(LUT_x,besseli(-1/4,LUT_x,1),'spline');
%         LUT_B34  = griddedInterpolant(LUT_x,besseli(3/4,LUT_x,1),'spline');
%         LUT_Bm34 = griddedInterpolant(LUT_x,besseli(-3/4,LUT_x,1),'spline');
%     elseif strcmp(SAMOSAimpl,'DPM_v2_5_2')
%         %Load look-up-tables alpha_p parameter and the F0 and F1 terms
%         fid      = fopen('LUT_Alpha_P.txt'); LUT_AP = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
%         LUT_AP   = griddedInterpolant(LUT_AP(:,1),LUT_AP(:,2),'linear','nearest');
%         fid      = fopen('F0.txt'); LUT_F0 = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
%         LUT_F0   = griddedInterpolant(LUT_F0(:,1),LUT_F0(:,2),'spline','nearest');
%         fid      = fopen('F1.txt'); LUT_F1 = cell2mat(textscan(fid,'%f %f\n','HeaderLines',1)); fclose(fid);
%         LUT_F1   = griddedInterpolant(LUT_F1(:,1),LUT_F1(:,2),'spline','nearest');
%     else
%         error('Which SAMOSA implementation should be used?')
%     end
% end
% 
% %Use DEM to compute initial retracking point
% Xpk_ALL = nan(size(CS.GEO.H(:)));
% if UseDEM2CompInitRP
%     Xpk_ALL = ((CS.GEO.H(:) - DEM(CS.GEO.LON(:),CS.GEO.LAT(:))) - CorrST01(:) + DDAcf.Biases.Range - CS.MEA.ref_range(:) + (0.5*DDAcf.c*RefBin*dt))/(0.5*DDAcf.c*dt);
%     Xpk_ALL(Xpk_ALL < 1 | Xpk_ALL > NrBins) = NaN;
% end
% 
% % profile on -detail builtin -history
% 
% %% Retrack waveforms
% for i = IDXpoi
%     %Copy normalized waveform i to vector WD
%     WD        = NORMdata(:,i);
%     
%     if ~IDX(i) || any(isnan(WD)); continue; end
%     if ~any(CS.WFc(i) == [0 1 4]); CS.HEI(i) = 100; continue; end
% 
% 
%     %Find largest peak in waveform and return associated bin index
%     [Ypk,Xpk] = max(WD(MINbin*os_ZP:MAXbin*os_ZP));
%     Xpk       = Xpk+(MINbin*os_ZP)-1;
%     if ~isnan(Xpk_ALL(i)), Xpk = Xpk_ALL(i); end
%     Wpk       = 10; %Just a number
% 
%     %Set initial values for SAMOSA retracker
%     % IDXrm     = max([1 i-10]):min([numel(CS.GEO.LAT) i+9]);
%     % [~,Xpk]   = max(prod(NORMdata(MINbin*os_ZP:MAXbin*os_ZP,IDXrm),2));
%     % Xpk       = Xpk+(MINbin*os_ZP)-1;
%     t0_0      = (Xpk-1 - RefBin)*dt*1E9;
% %     IDXrm     = max([1 i-20]):max([1 i-1]);
% %     Pu0       = nanmean(CS.Pu(IDXrm));
% %     SWH0      = nanmean(CS.SWH(IDXrm));
% %     nu0       = nanmean(CS.nu(IDXrm));
%     Pu0       = 1;
%     SWH0      = 2;
%     nu0       = 10;
%     
%     %Apply retracking
%     switch Retracker
%         case 'BetaX'
%             %X-parameter Beta-retracker with a linear trailing edge (Martin
%             %et al., 1983)
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'linear');
%         case 'BetaX_expTE'
%             %X-parameter Beta-retracker with an exponential trailing edge
%             %(Deng & Featherstone, 2006)
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'exponential');
%         case 'Brown'
%             %Brown Theoretical Ocean Model (Passaro et al., 2014)
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Brown(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
%         case 'BrownHayne'
%             %Brown-Hayne Theoretical Ocean Model (Gommenginger et al., 2011)
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BrownHayne(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
%         case 'D2P'
%             %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
%             %(Giles et al., 2007)
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.D2P(WD,BinIDs,Ypk,Xpk,Wpk);
%         case 'FunctionFit'
%             %"Function Fit" retracker (Surface Topography Mission (STM)
%             %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
%             %[SD-03] [SD-07])
%             [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.FunctionFit(WD,BinIDs,Ypk,Xpk,Wpk);
%         case 'OCOG'
%             %Offset Centre Of Gravity retracker (Wingham et al., 1986)
%             [CS.n_ret(i),CS.Pu(i)] = RETRACKER.OCOG(WD,BinIDs');
%         case 'Threshold'
%             %Threshold retracker (Davis, 1997).
%             [CS.n_ret(i),CS.Pu(i)] = RETRACKER.Threshold(WD,BinIDs',WD);
%         case 'SAMOSA2'
%             %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
%             [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,CS.SAR.BeamAngle{i}');
%         case 'SAMOSA2FF'
%             %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
%             [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,0.5*pi);
%         case 'SAMOSA3'
%             %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
%             [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,CS.SAR.BeamAngle{i}');
%         case 'SAMOSA3FF'
%             %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
%             [x0,lb,ub,XDATA,IDXnr,BeamIDX] = RETRACKER.SAMOSA(WD,BinIDs,DDAcf,CS.MEA.TN(i),Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),CS.GEO.H(i),CS.GEO.V.V(i),CS.GEO.orbit_slope(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,0.5*pi);
%         otherwise
%             error('Retracker not implemented!')
%     end
% 
%     %Set function handle
%     if strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'SD')
%         fun = @(x,n) RETRACKER.SAMOSAfun(x,n,DDAcf,LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34,BeamIDX);
%     elseif strncmp(Retracker,'SAMOSA',6) && strcmp(SAMOSAimpl,'DPM_v2_5_2')
%         fun = @(x,n) RETRACKER.SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX);
%     end
%     
%     %In case of analytical retrackers, solve the non-linear obs. eq.
%     if ~any(strcmp(Retracker,{'OCOG','Threshold'}))
%         %In case of the SAMOSA retracker, the 4 unknown parameters are not
%         %solved simultaneously. For ocean waveforms, nu (the inverse of the
%         %mean-square slope of the sea surface) is set to 0 and not
%         %estimated. In case of lead waveforms, the SWH is set to 0 and not
%         %estimated. For land contaminated waveforms, Dinardo et al. (2018)
%         %apply a dual step retracking where first Pu, t0, and SWH are
%         %estimated. Thereafter, the SWH is set 0 and Pu, t0, and nu are
%         %(re-)estimated.
%         
%         %Solve for parameters in case of ocean waveform
%         if any(CS.WFc(i) == [0 1 5])
%             %Estimate Pu, t0, and SWH
%             XDATA(end)           = 3;
%             x0_tmp               = x0(1:3); lb_tmp = lb(1:3); ub_tmp = ub(1:3);
%             if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
%             [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
%             x(4)                 = 0;
% 
%             %Verify reason why algorithm is terminated
%             CS.ExitF(i)  = exitflag;
%             if exitflag <= 0, continue, end
%             
%             %Assess whether waveform is land contaminated based on quality
%             %of fit
%             CS.PCorr(i) = corr(WD,fun(x,XDATA),'type','Pearson');
%             if CS.PCorr(i) <= 0.9, CS.WFc(i) = 6; XDATA(end-8) = x(3); end
%             %if 100*sqrt(Resnorm/NrBins) > 4, CS.WFc(i) = 6; XDATA(end-8) = x(3); end
%         end
% 
%         %Solve for parameters in case of lead or land contaminated waveform
%         if any(CS.WFc(i) == [4 6])
%             if CS.WFc(i) == 4, XDATA(end-8) = 0; end
%             
%             %(Re-)estimate Pu, t0, and nu (inverse of the mean-square slope of
%             %the sea surface).
%             XDATA(end)           = 4;
%             x0_tmp               = x0([1 2 4]); lb_tmp = lb([1 2 4]); ub_tmp = ub([1 2 4]);
%             if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
%             [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
%             x                    = [x(1:2) XDATA(end-8) x(3)];
% 
%             %Verify again the reason why algorithm is terminated
%             CS.ExitF(i) = exitflag;
%             if exitflag <= 0, continue, end
% 
%             %Assess whether waveform has indeed a specular shape
%             CS.PCorr(i) = corr(WD,fun(x,XDATA),'type','Pearson');
%             if CS.PCorr(i) <= 0.95, CS.WFc(i) = 99; CS.PCorr(i) = NaN; end
%         end
% 
%         %Apply threshold retracker in case waveform belongs to classes 2,
%         %3, or 99
%         if any(CS.WFc(i) == [2 3 99])
%             x           = nan(1,4);
%             [x(2),x(1)] = RETRACKER.Threshold(WD,BinIDs',WD);
%             x(2)        = (x(2) - RefBin)*(dt*1E9);
%             Resnorm     = NaN;
%         end
%         
%         %Copy retracking location [bins] and other estimated parameters to CS
%         if any(strcmp(Retracker,{'Brown','BrownHayne'}))
%             CS.n_ret(i) = x(IDXnr)/(dt*1E9);
%             CS.Pu(i)    = x(IDXnr-1);
%         elseif any(strcmp(Retracker,{'SAMOSA2','SAMOSA2FF','SAMOSA3','SAMOSA3FF'}))
%             CS.n_ret(i) = (x(IDXnr)/(dt*1E9)) + RefBin;
%             CS.Pu(i)    = x(IDXnr-1);
%             CS.SWH(i)   = x(IDXnr+1);
%             CS.nu(i)    = x(IDXnr+2);
%             CS.MF(i)    = 100*sqrt(Resnorm/NrBins);
%         else
%             CS.n_ret(i) = x(IDXnr);
%             CS.Pu(i)    = x(IDXnr-1);
%         end
%         
%         %Reconstruct waveform
%         CS.WDrecon(:,i) = fun(x,XDATA);
%     end
% end
% clear('WD','sub_n','sub_WD','Ypk','Xpk','IDXpks','i','j')
% 
% % profile viewer
% % profile off
% 
% %% Data editing
% CS.n_ret(CS.n_ret < MINbin*os_ZP | CS.n_ret > MAXbin*os_ZP) = NaN;
% 
% %% Compute corrected range (corrections include instrumental, propagation and geophysical corrections)
% %Compute range Eqn 2.8�?1 CryoSat Product Handbook
% CS.('range')      = CS.MEA.ref_range + (0.5*DDAcf.c*(CS.n_ret*dt - RefBin*dt)) - DDAcf.Biases.Range;
% CS.('range_SSHi') = CS.range;
% 
% %Apply propagation/geophysical corrections
% ST                     = CS.surf_type;
% CS.('SumCorrST')       = nan(size(CorrST01)); CS.('FlgC')       = true(size(CorrST01));
% CS.SumCorrST(ST <= 1)  = CorrST01(ST <= 1);   CS.FlgC(ST <= 1)  = FlgC01(ST <= 1);
% CS.SumCorrST(ST >= 2)  = CorrST23(ST >= 2);   CS.FlgC(ST >= 2)  = FlgC23(ST >= 2);
% CS.('SumCorrSTi')      = CS.SumCorrST;        CS.('FlgCi')      = CS.FlgC;
% CS.SumCorrSTi(ST <= 1) = CorrST01i(ST <= 1);  CS.FlgCi(ST <= 1) = FlgC01i(ST <= 1);
% CS.range               = CS.range + CS.SumCorrST;
% CS.range_SSHi          = CS.range_SSHi + CS.SumCorrSTi;
% clear('CorrST01','CorrST01i','CorrST23','ST','FlgC01','FlgC01i','FlgC23')
% 
% %Compute (instantaneous) (sea) surface heights relative to the reference ellipsoid
% CS.('HEI')  = CS.GEO.H(:) - CS.range;
% CS.('SSHi') = CS.GEO.H(:) - CS.range_SSHi;
% 
% %% Compute backscatter coefficient (sigma0 / sigma-naught)
% switch SAT
%     case {'CS','CryoSat'}
%         CS.Pu         = CS.Pu ./ NORMfactor(:);
%         CS.('sigma0') = Compute_Sigma0(CS,DDAcf);
%         CS.Pu         = CS.Pu .* NORMfactor(:);
%     case {'S3A','Sentinel-3A','S3B','Sentinel-3B'}
%         %Surface Topography Mission (STM) SRAL/MWR L2 Algorithms
%         %Definition, Accuracy and Specification, Section 2.15.3.3
%         CS.('sigma0') = CS.SAR.scale_factor_ku_l1b_echo_sar_ku + 10*log10(CS.Pu) + DDAcf.Biases.sigma0;
%     otherwise
%         error('SAT: %s not recognized',SAT)
% end
% 
% %% Analyze output
% if MakePlots
%     for i = IDXpoi
%         if all(isnan(CS.n_ret(i,:))), continue, end
%     
%         figure('Position',get(0,'Screensize'));
% 
%         %Select valid entries
%         IDXvalid = ~isnan(CS.n_ret(i,:));
% 
%         %Plot observed/reconstructed waveform
%         WD        = CS.SAR.data(:,IDXpoi);
%         WDnorm    = NORMdata(:,i);
%         StrYlabel = 'Power (Watts)';
%         if ~any(strcmp(Retracker,{'OCOG','Threshold'})); WD = WDnorm; StrYlabel = 'Normalized power'; end
%         subplot(2,3,1),plot(1:NrBins,WD,'.-','LineWidth',2),hold on
%         if ~all(isnan(CS.WDrecon(:,IDXpoi)))
%             plot(1:NrBins,CS.WDrecon(:,IDXpoi),'g-','LineWidth',2)
%             plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,CS.WDrecon(:,IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
%             LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
%         else
%             plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
%             LegendEntries = {'Observed waveform','Retracking points'};
%         end
%         axis([0 NrBins 0 max(WD)])
%         title('Observed waveform','fontsize',FntSz)
%         legend(LegendEntries,'box','off')
%         xlabel('Bin','fontsize',FntSz)
%         ylabel(StrYlabel,'fontsize',FntSz)
%         set(gca,'fontsize',FntSz)
% 
%         WD        = CS.SAR.data(:,IDXpoi);
%         subplot(2,3,2),plot(1:NrBins,WD,'.-','LineWidth',2),hold on
%         if ~all(isnan(CS.WDrecon(:,IDXpoi)))
%             plot(1:NrBins,fun(x,XDATA)/NORMfactor(IDXpoi),'g-','LineWidth',2)
%             plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,fun(x,XDATA)/NORMfactor(IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
%             LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
%         else
%             plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
%             LegendEntries = {'Observed waveform','Retracking points'};
%         end
%         plot(1:NrBins,ones(1,NrBins)*CS.Pu(IDXpoi)/NORMfactor(IDXpoi),'k--');
%         LegendEntries = horzcat(LegendEntries,{'Pu'});
%         axis([0 NrBins 0 max(WD)])
%         title('Observed waveform','fontsize',FntSz)
%         legend(LegendEntries,'box','off')
%         xlabel('Bin','fontsize',FntSz)
%         ylabel('Power (Watts)','fontsize',FntSz)
%         set(gca,'fontsize',FntSz)
% 
%         %Plot coastline
%         LOLA = Get_River_Border_data('GSHHS','coastlines','f',DOM(1,:),DOM(2,:));
%         subplot(2,3,4),plot(LOLA(:,1),LOLA(:,2),'k-'),hold on
%         plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.')
%         plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
%         axis([minmax(LOLA(:,1)) minmax(LOLA(:,2))])
%         xlabel('Longitude (deg)','fontsize',FntSz)
%         ylabel('Latitude (deg)','fontsize',FntSz)
%         set(gca,'fontsize',FntSz)
% 
%         %Plot DEM around point of interest including ground track
%         if ~isempty(DEM)
%         AZIM   = azimuth(CS.GEO.LAT(1:end-1),CS.GEO.LON(1:end-1),CS.GEO.LAT(2:end),CS.GEO.LON(2:end),DDAcf.RefEll);
%         MaxDoN = km2deg(20,6371); %Maximum distance-of-nadir [km]
%         DEMx   = DEM.GridVectors{1}; DEMy = DEM.GridVectors{2};
%         [~,Ix] = min(abs(DEMx-CS.GEO.LON(i))); [~,Iy] = min(abs(DEMy-CS.GEO.LAT(i))); 
%         MDidx  = ceil(MaxDoN/mode(diff(DEMx)));
%         subplot(2,3,5),imagesc(DEMx(Ix+(-MDidx:MDidx)),DEMy(Iy+(-MDidx:MDidx)),DEM.Values((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx)))),hold on
%         axis square, colorbar, set(gca,'YDir','normal')
%         set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String','meters','FontSize',FntSz);
%         %Plot coastline
%         plot(LOLA(:,1),LOLA(:,2),'-','Color',[.5 .5 .5]),hold on
%         %Plot direction of flight
%         quiver(CS.GEO.LON(i),CS.GEO.LAT(i),(CS.GEO.LON(i+1)-CS.GEO.LON(i)),(CS.GEO.LAT(i+1)-CS.GEO.LAT(i)),20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
%         %Plot across-track points
%         [lat_ac,lon_ac] = reckon(CS.GEO.LAT(i),CS.GEO.LON(i),deg2km(-MaxDoN:mode(diff(DEMx))/2:MaxDoN,6371)*1000,AZIM(i)-90,DDAcf.RefEll);
%         plot(lon_ac,lat_ac,'k.-')
%         %Plot identified scatterers
%         plot(CS.GEO.LON(i),CS.GEO.LAT(i),'o','Color',[.5 .5 .5],'LineWidth',2,'MarkerSize',12)
%         title('DEM','fontsize',FntSz)
%         xlabel('Longitude (deg)','fontsize',FntSz)
%         ylabel('Latitude (deg)','fontsize',FntSz)
%         set(gca,'fontsize',FntSz)
%         end
% 
%         % subplot(2,3,6),plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.'),hold on
%         % plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
%         % axis([minmax(DEMx(Ix+(-MDidx:MDidx))) minmax(DEMy(Iy+(-MDidx:MDidx)))])
%         % plot_google_map
%         % xlabel('Longitude (deg)','fontsize',FntSz)
%         % ylabel('Latitude (deg)','fontsize',FntSz)
%         % set(gca,'fontsize',FntSz)
% 
%         %Set background color to white
%         set(gcf, 'Color',[1 1 1])
%     end
% end

%% Save output
DATA              = struct;
%IDX               = ~isnan(CS.HEI);
DATA.('TIME')     = CS.TIME(IDX);
DATA.('LAT')      = CS.GEO.LAT(IDX);
DATA.('LON')      = CS.GEO.LON(IDX);
%DATA.('HEI')      = single(CS.HEI(IDX));
%DATA.('sumCORR')  = single(CS.SumCorrST(IDX));
% DATA.('SAL')      = CS.COR.SAL(IDX);
%DATA.('SSB')      = 
%DATA.('FlgC')     = CS.FlgC(IDX);
%DATA.('SSHi')     = single(CS.SSHi(IDX));
%DATA.('sumCORRi') = single(CS.SumCorrSTi(IDX));
%DATA.('FlgCi')    = CS.FlgCi(IDX);
DATA.('sigma0')   = single(CS.sigma0(IDX));
%DATA.('SWH')      = single(CS.SWH(IDX));
%DATA.('SurfT')    = int8(CS.surf_type(IDX));
%DATA.('WFc')      = int8(CS.WFc(IDX));
%DATA.('nu')       = single(CS.nu(IDX));
%DATA.('ExitF')    = int8(CS.ExitF(IDX));
%DATA.('MF')       = single(CS.MF(IDX));
%DATA.('PCorr')    = single(CS.PCorr(IDX));
%DATA.('cEqTide')  = single(cEqTide(IDX));
%DATA.('cLPTide')  = single(cLPTide(IDX));
%DATA.('cOLTide')  = single(cOLTide(IDX));
%DATA.('cSETide')  = single(cSETide(IDX));
%DATA.('cGPTide')  = single(cGPTide(IDX));
% DATA.('cHF')      = single(cHF(IDX));
%DATA.('cMeteo')  = single(cMeteo(IDX));

end
