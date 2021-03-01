function [DATA,CS] = CryoSat_SAR_L1b_to_L2(FName,DOM,Retracker,SetRetr,DEM,IDXpoi)

%CRYOSAT_SAR_L1B_TO_L2 processes level 2 data from CryoSat baseline B/C
%level 1b SAR data.

%% Settings
%General settings
LoadCommonSettings

%Computation settings
defval('SetRetr',{})
for i=1:size(SetRetr,1), eval(sprintf('%s = %f;',SetRetr{i,1},SetRetr{i,2})); end
defval('FName','2016/03/CS_OFFL_SIR_SAR_1B_20160302T183427_20160302T183612_C001') %*.DBL file that contains level 1b data
defval('DOM',[59 84;-74 -10])                                                     %String that identifies region of interest/geographic quadrangle specified by [minLAT maxLAT; minLON maxLON]
defval('Retracker','Threshold')                                                   %Retracker to be used
defval('SolverAlgo','TRR')                                                        %Solver Algorithm to be applied in case Retracker is an analytical retracker ('LM' = Levenberg-Marquardt; 'TRR' = trust-region-reflective)
defval('max_waveform_energy',150)                                                 %Maximum allowed energy of NORMALIZED waveform (i.e., area below the curve)
defval('MAfromSTR',true)                                                          %Obtain mispointing angle from star tracker data

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

%Physical constants & SIRAL Characteristics
c       = 299792458;                       %Light velocity [m/s]
B       = 320E6;                           %Pulse bandwidth [Hz]
WGS84   = referenceEllipsoid('wgs84','m'); %Reference ellipsoid used for CryoSat

%Set paths
PathL1b = fullfile(PathDATA,'RadAlt','CryoSat','SIR_SAR_L1'); %Path to data

%Remaining settings
defval('MakePlots',exist('IDXpoi','var')) %Make plots
if MakePlots; if ~isscalar(IDXpoi); MakePlots = false; end; end
defval('FntSz',14)                        %Set fontsize to be used in figures

%% Read & Crop CryoSat Level 1b BASELINE C SAR data
%Read data
if isstruct(FName)
    CS     = FName;
else
    [~,CS] = Cryo_L1b_read(fullfile(PathL1b,sprintf('%s.DBL',FName)));
end

%Remove field AVG to save RAM
CS = rmfield(CS,'AVG');

%Correct window delay for DORIS USO drift (pp 20 CryoSat Product Handbook)
CS.MEA.win_delay = CS.MEA.win_delay.*(CS.GEO.USO+1);

%Set baseline dependent constants to be used in processing
%https://wiki.services.eoportal.org/tiki-index.php?page=CryoSat+Technical+Notes
% ARESYS - Known biases in CryoSat Level1b products
% ARESYS - Main evolutions and expected quality improvements in Baseline C Level 1b products
if strcmp(CS.GEO.Baseline_ID,'C')
    NrBins = 256;             %Nr of bins/samples in any waveform
    dt     = 1/(2*B);         %Waveform sampling interval [s]
    defval('MINbin',20);      %Minimum bin index of interval in which peaks are detected
    defval('MAXbin',250);     %Maximum bin index of interval in which peaks are detected

    %Residual errors        = [applied, residual]
    BaselineInfo.Datation   = [-0.5195,-0.029];  % [ms]
    BaselineInfo.Range      = [0.6730,  47E-3];  % [m]
else
    error('Baseline_ID not valid')
end

%Get Roll, Yaw, and Pitch angles from star tracker data
if MAfromSTR
    try
        tSTR             = datetime(CS.GEO.Start_Time./(24.*60.*60) + datenum('01-Jan-2000 00:00:00') - 1,'ConvertFrom','datenum');
        mpSTR            = READ_STR_Mispointing_Angles(tSTR);
        tOBS             = datenum('2000','yyyy') + CS.GEO.TAI.days + CS.GEO.TAI.secs./86400 + CS.GEO.TAI.microsecs./1e6./86400;
        DUM              = interp1(mpSTR.Time,mpSTR.Roll,tOBS,'spline');
        IDXnan           = isnan(interp1(mpSTR.Time,mpSTR.Roll,tOBS));
        DUM(IDXnan)      = NaN;
        CS.GEO.Antenna_Bench_Roll(~isnan(DUM)) = DUM(~isnan(DUM));
        DUM              = interp1(mpSTR.Time,mpSTR.Yaw,tOBS,'spline');
        IDXnan           = isnan(interp1(mpSTR.Time,mpSTR.Yaw,tOBS));
        DUM(IDXnan)      = NaN;
        %Note that yaw angle is defined in different sign convention!
        CS.GEO.Antenna_Bench_Yaw(~isnan(DUM)) = -DUM(~isnan(DUM));
        DUM              = interp1(mpSTR.Time,mpSTR.Pitch,tOBS,'spline');
        IDXnan           = isnan(interp1(mpSTR.Time,mpSTR.Pitch,tOBS));
        DUM(IDXnan)      = NaN;
        %Note that pitch angle is defined in different sign convention!
        CS.GEO.Antenna_Bench_Pitch(~isnan(DUM)) = -DUM(~isnan(DUM));
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

% %Compute mispointing angles
% CS.GEO.('MPA') = acos(cos(CS.GEO.Antenna_Bench_Pitch)./sqrt(sum(cat(3,sin(CS.GEO.Antenna_Bench_Roll),sin(CS.GEO.Antenna_Bench_Yaw),cos(CS.GEO.Antenna_Bench_Pitch)).^2,3)));

%Compute the local radii of curvature of the Earth's surface (Maulik Jain,
%Improved sea level determination in the Arctic regions through development
%of tolerant altimetry retracking, Eq. 5.2, pp. 47)
CS.GEO.('Re')  = sqrt(WGS84.SemimajorAxis^2*cosd(CS.GEO.LAT).^2 + WGS84.SemiminorAxis^2*sind(CS.GEO.LAT).^2);

%The power echo	sample values are all scaled to	fit	between	0 and 65535.
%The scaling factors can change	for	each waveform. To convert these back to
%values in Watts the following equation should be used (CryoSat Product
%Handbook, Eqn 4.2‐1):
%Power in Watts	= scaled value * (scale	factor * 10^‐9) * 2^scale power
CS.SAR.data = bsxfun(@times,CS.SAR.data,reshape((CS.SAR.echo_scaling.*1E-9) .* 2.^(CS.SAR.echo_scale_power),1,size(CS.SAR.data,2),size(CS.SAR.data,3)));

%Normalize each waveform by max value of waveform
NORMfactor = 1./max(CS.SAR.data,[],1);
NORMdata   = NORMfactor .* CS.SAR.data;

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
%Identify flagged data
IDXfd = CS.GEO.MCD_FLAG.Block_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Blank_Block(:) == 1 | ...
    CS.GEO.MCD_FLAG.Datation_Degraded(:) == 1 | CS.GEO.MCD_FLAG.Orbit_Propag_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.Echo_Saturation(:) == 1 | CS.GEO.MCD_FLAG.Other_Echo_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.Rx1_Err_SARin(:) == 1 | CS.GEO.MCD_FLAG.Rx2_Err_SARin(:) == 1 | ...
    CS.GEO.MCD_FLAG.Wind_Delay_Incon(:) == 1 | CS.GEO.MCD_FLAG.AGC_Incon(:) == 1 | ...
    CS.GEO.MCD_FLAG.TRK_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.RX1_ECHO_Err(:) == 1 | ...
    CS.GEO.MCD_FLAG.RX2_ECHO_Err(:) == 1 | CS.GEO.MCD_FLAG.NPM_Incon(:) == 1 | ...
    CS.GEO.MCD_FLAG.Power_Scaling_Err(:) == 1;
IDX(IDXfd) = false;

%Identify waveforms for which power == 0 for all entries
IDX(squeeze(all(CS.SAR.data == 0,1))) = false;

%Energy NORMALIZED waveform < max_waveform_energy && power at first bins of
%NORMALIZED waveform should be at noise level
IDX(squeeze(trapz(1:NrBins,NORMdata,1) >= max_waveform_energy) | squeeze(any(NORMdata(1:MINbin,:,:) > .1))) = false;

%Return if no data remain
if ~any(IDX)
    DATA = struct('TIME',[],'LAT',[],'LON',[],'HEI',[],'SurfT',[]);
    return
end

%Select points of interest (by default entire track is retracked)
defval('IDXpoi',find(IDX)')

%% Classify waveforms
%Preliminaries
[CS.('n_ret'),CS.('Pu'),CS.('WFc')] = deal(nan(numel(CS.GEO.LAT),1));

%Apply threshold retracker to obtain retracking point and Pu
[CS.n_ret(IDXpoi),CS.Pu(IDXpoi)]    = RETRACKER.Threshold_mat(NORMdata(:,IDXpoi));

%Compute range and backscatter coefficient (sigma0 / sigma-naught)
CS.('range')   = (0.5*c*CS.MEA.win_delay(:)) + (0.5*c*(CS.n_ret*dt - ((NrBins/2)+1)*dt)) - BaselineInfo.Range(2);
CS.Pu          = CS.Pu ./ NORMfactor(:);
CS.ClassPar.sigma0 = Compute_Sigma0(CS)';
CS.Pu          = CS.Pu .* NORMfactor(:);

%Classify waveforms
[CS.WFc(IDXpoi),CS.ClassPar.Ptail(IDXpoi), CS.ClassPar.PP(IDXpoi),CS.ClassPar.RatioIP(IDXpoi),CS.ClassPar.kurt(IDXpoi),CS.ClassPar.PPloc(IDXpoi)] = Classify_Waveforms(NORMdata(:,IDXpoi),CS.ClassPar.sigma0(IDXpoi));

%% Retrack waveforms
%Preliminaries
[CS.('n_ret'),CS.('Pu'),CS.('SWH'),CS.('sigma0')] = deal(nan(numel(CS.GEO.LAT),1));
[CS.('nu'),CS.('ExitF'),CS.('MF')]  = deal(nan(numel(CS.GEO.LAT),1));
BinIDs                              = (1:NrBins)';
CS.('WDrecon')                      = nan(size(CS.SAR.data));

%Generate look-up-tables for fast evaluation of modified Bessel functions
%of the first kind
LUT_x    = logspace(-16,4,1000)';
LUT_B14  = griddedInterpolant(LUT_x,besseli(1/4,LUT_x,1),'spline');
LUT_Bm14 = griddedInterpolant(LUT_x,besseli(-1/4,LUT_x,1),'spline');
LUT_B34  = griddedInterpolant(LUT_x,besseli(3/4,LUT_x,1),'spline');
LUT_Bm34 = griddedInterpolant(LUT_x,besseli(-3/4,LUT_x,1),'spline');

% profile on -detail builtin -history

%Apply retracking
for i = IDXpoi(:)'
    if ~IDX(i); continue; end
    if CS.WFc(i) ~= 2; continue; end
    
    %Copy normalized waveform i to vector WD
    WD        = NORMdata(:,i);

    %Find largest peak in waveform and return associated bin index
    [Ypk,Xpk] = max(WD(MINbin:MAXbin));
    Xpk       = Xpk+MINbin-1;
    Wpk       = 10; %Just a number

    %Set initial values for SAMOSA retracker
    % IDXrm     = max([1 i-10]):min([numel(CS.GEO.LAT) i+9]);
    % [~,Xpk]   = max(prod(NORMdata(MINbin:MAXbin,IDXrm),2));
    % Xpk       = Xpk+MINbin-1;
    t0_0      = (Xpk-1 - ((NrBins/2)+1))*dt*1E9;
%     IDXrm     = max([1 i-20]):max([1 i-1]);
%     Pu0       = nanmean(CS.Pu(IDXrm));
%     SWH0      = nanmean(CS.SWH(IDXrm));
%     nu0       = nanmean(CS.nu(IDXrm));
    Pu0       = 1;
    SWH0      = 2;
    nu0       = 10;
    
    %Apply retracking
    switch Retracker
        case 'BetaX'
            %X-parameter Beta-retracker with a linear trailing edge (Martin
            %et al., 1983)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'linear');
        case 'BetaX_expTE'
            %X-parameter Beta-retracker with an exponential trailing edge
            %(Deng & Featherstone, 2006)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Beta5(WD,BinIDs,Ypk,Xpk,'exponential');
        case 'Brown'
            %Brown Theoretical Ocean Model (Passaro et al., 2014)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.Brown(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'BrownHayne'
            %Brown-Hayne Theoretical Ocean Model (Gommenginger et al., 2011)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.BrownHayne(WD,BinIDs,Ypk,Xpk,CS.GEO.H(i),dt*1E9);
        case 'D2P'
            %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
            %(Giles et al., 2007)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.D2P(WD,BinIDs,Ypk,Xpk,Wpk);
        case 'FunctionFit'
            %"Function Fit" retracker (Surface Topography Mission (STM)
            %SRAL/MWR L2 Algorithms Definition, Accuracy and Specification
            %[SD-03] [SD-07])
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.FunctionFit(WD,BinIDs,Ypk,Xpk,Wpk);
        case 'OCOG'
            %Offset Centre Of Gravity retracker (Wingham et al., 1986)
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.OCOG(WD,BinIDs');
        case 'Threshold'
            %Threshold retracker (Davis, 1997).
            [CS.n_ret(i),CS.Pu(i)] = RETRACKER.Threshold(WD,BinIDs',WD);
        case 'SAMOSA2'
            %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.SAMOSA(WD,BinIDs,Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),0.5*c*CS.MEA.win_delay(i),CS.GEO.V.V(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,floor(CS.SAR.N_averaged_echoes(i)/4),LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34);
        case 'SAMOSA2FF'
            %SAMOSA2 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.SAMOSA(WD,BinIDs,Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),0.5*c*CS.MEA.win_delay(i),CS.GEO.V.V(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),1,0,1,LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34);
        case 'SAMOSA3'
            %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.SAMOSA(WD,BinIDs,Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),0.5*c*CS.MEA.win_delay(i),CS.GEO.V.V(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,floor(CS.SAR.N_averaged_echoes(i)/4),LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34);
        case 'SAMOSA3FF'
            %SAMOSA3 retracker (Ray et al. 2015; Dinardo et al., 2018)
            [fun,x0,lb,ub,XDATA,IDXnr] = RETRACKER.SAMOSA(WD,BinIDs,Pu0,t0_0,SWH0,nu0,CS.GEO.Re(i),0.5*c*CS.MEA.win_delay(i),CS.GEO.V.V(i),CS.GEO.Antenna_Bench_Pitch(i),CS.GEO.Antenna_Bench_Roll(i),0,0,1,LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34);
        otherwise
            error('Retracker not implemented!')
    end

    %In case of analytical retrackers, solve the non-linear obs. eq.
    if ~any(strcmp(Retracker,{'OCOG','Threshold'}))
        %In case of the SAMOSA retracker, the 4 unknown parameters are not
        %solved simultaneously. For ocean waveforms, nu (the inverse of the
        %mean-square slope of the sea surface) is set to 0 and not
        %estimated. In case of lead waveforms, the SWH is set to 0 and not
        %estimated. For land contaminated waveforms, Dinardo et al. (2018)
        %apply a dual step retracking where first Pu, t0, and SWH are
        %estimated. Thereafter, the SWH is set 0 and Pu, t0, and nu are
        %(re-)estimated.
        
        %Solve for parameters in case of ocean waveform
        if any(CS.WFc(i) == [0 1 3])
            %Estimate Pu, t0, and SWH
            XDATA(end-1)         = 3;
            x0_tmp               = x0(1:3); lb_tmp = lb(1:3); ub_tmp = ub(1:3);
            if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
            [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
            x(4)                 = 0;

            %Verify reason why algorithm is terminated
            CS.ExitF(i)  = exitflag;
            if exitflag <= 0, continue, end
            
            %Assess whether or not ocean waveform is a land contaminated
            %waveform
            if 100*sqrt(Resnorm/NrBins) > 4, CS.WFc(i) = 4; XDATA(end-8) = x(3); end
        end

        %Solve for parameters in case of lead or land contaminated waveform
        if any(CS.WFc(i) == [2 4])
            if CS.WFc(i) == 2, XDATA(end-8) = 0; end
            
            %(Re-)estimate Pu, t0, and nu (inverse of the mean-square slope of
            %the sea surface).
            XDATA(end-1)         = 4;
            x0_tmp               = x0([1 2 4]); lb_tmp = lb([1 2 4]); ub_tmp = ub([1 2 4]);
            if strcmp(SolverAlgo,'LM'), lb_tmp = []; ub_tmp = []; end
            [x,Resnorm,exitflag] = DARTutils.SolveNLLS(fun,x0_tmp,XDATA,WD,lb_tmp,ub_tmp,options);
            x                    = [x(1:2) XDATA(end-8) x(3)];

            %Verify again the reason why algorithm is terminated
            CS.ExitF(i) = exitflag;
            if exitflag <= 0, continue, end
        end
        
        %Copy retracking location [bins] and other estimated parameters to CS
        if any(strcmp(Retracker,{'Brown','BrownHayne'}))
            CS.n_ret(i) = x(IDXnr)/(dt*1E9);
            CS.Pu(i)    = x(IDXnr-1);
        elseif any(strcmp(Retracker,{'SAMOSA2','SAMOSA2FF','SAMOSA3','SAMOSA3FF'}))
            CS.n_ret(i) = (x(IDXnr)/(dt*1E9)) + ((NrBins/2)+1);
            CS.Pu(i)    = x(IDXnr-1);
            CS.SWH(i)   = x(IDXnr+1);
            CS.nu(i)    = x(IDXnr+2);
            CS.MF(i)    = 100*sqrt(Resnorm/NrBins);
        else
            CS.n_ret(i) = x(IDXnr);
            CS.Pu(i)    = x(IDXnr-1);
        end
        
        %Reconstruct waveform
        CS.WDrecon(:,i) = fun(x,XDATA);
    end
end
clear('WD','sub_n','sub_WD','Ypk','Xpk','IDXpks','i','j')

% profile viewer
% profile off

%% Data editing
CS.n_ret(CS.n_ret < MINbin | CS.n_ret > MAXbin) = NaN;

%% Compute corrected range (range corrections include instrumental, range and geophysical corrections)
%Compute range Eqn 2.8‐1 CryoSat Product Handbook
CS.('range')      = (0.5*c*CS.MEA.win_delay(:)) + (0.5*c*(CS.n_ret*dt - ((NrBins/2)+1)*dt)) - BaselineInfo.Range(2);
CS.('range_SSHi') = CS.range;
CS.('rangebe4corr') = CS.range;

%Corrections to be applied in case (i) surface == open oceans or
%semi‐enclosed seas OR (ii) surface == enclosed seas or lakes (SSB corr. is
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

%Apply corrections by nearest neighbour interpolation
ST                     = reshape(repmat(CS.COR.surf_type,20,1),numel(CS.COR.surf_type)*20,1);
CorrST01               = reshape(repmat(CorrST01,20,1),numel(CorrST01)*20,1);
CorrST01i              = reshape(repmat(CorrST01i,20,1),numel(CorrST01i)*20,1);
CorrST23               = reshape(repmat(CorrST23,20,1),numel(CorrST23)*20,1);
CS.range(ST <= 1)      = CS.range(ST <= 1) + CorrST01(ST <= 1);
CS.range(ST >= 2)      = CS.range(ST >= 2) + CorrST23(ST >= 2);
CS.range_SSHi(ST <= 1) = CS.range_SSHi(ST <= 1) + CorrST01i(ST <= 1);
CS.range_SSHi(ST >= 2) = CS.range_SSHi(ST >= 2) + CorrST23(ST >= 2);
CS.('surf_type')       = reshape(repmat(CS.COR.surf_type,20,1),numel(CS.COR.surf_type)*20,1);
CS.('CorrST01i')       = CorrST01i;
clear('CorrST01','CorrST01i','CorrST23','ST')

%Compute (instantaneous) (sea) surface heights relative to the WGS84 ellipsoid
CS.('HEI')  = CS.GEO.H(:) - CS.range;
CS.('SSHi') = CS.GEO.H(:) - CS.range_SSHi;

%Transform acquisition time to datenum format
CS.('TIME') = datenum('2000','yyyy') + CS.GEO.TAI.days(:) + CS.GEO.TAI.secs(:)./86400 + CS.GEO.TAI.microsecs(:)./1e6./86400 - BaselineInfo.Datation(2)/1E3/86400;

%% Compute backscatter coefficient (sigma0 / sigma-naught)
CS.Pu         = CS.Pu ./ NORMfactor(:);
CS.('sigma0') = Compute_Sigma0(CS);
CS.Pu         = CS.Pu .* NORMfactor(:);

%% Analyze output
if MakePlots
    for i = IDXpoi
        if all(isnan(CS.n_ret(i,:))), continue, end
    
        figure('Position',get(0,'Screensize'));

        %Select valid entries
        IDXvalid = ~isnan(CS.n_ret(i,:));

        %Plot observed/reconstructed waveform
        WD        = CS.SAR.data(:,IDXpoi);
        WDnorm    = NORMdata(:,i);
        StrYlabel = 'Power (Watts)';
        if ~any(strcmp(Retracker,{'OCOG','Threshold'})); WD = WDnorm; StrYlabel = 'Normalized power'; end
        subplot(2,3,1),plot(1:NrBins,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:NrBins,CS.WDrecon(:,IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,CS.WDrecon(:,IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        axis([0 NrBins 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel(StrYlabel,'fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        WD        = CS.SAR.data(:,IDXpoi);
        subplot(2,3,2),plot(1:NrBins,WD,'.-','LineWidth',2),hold on
        if ~all(isnan(CS.WDrecon(:,IDXpoi)))
            plot(1:NrBins,fun(x,XDATA)/NORMfactor(IDXpoi),'g-','LineWidth',2)
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,fun(x,XDATA)/NORMfactor(IDXpoi),CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Reconstructed waveform','Retracking points'};
        else
            plot(CS.n_ret(i,1:nnz(IDXvalid)),interp1(1:NrBins,WD,CS.n_ret(i,1:nnz(IDXvalid))),'ro','LineWidth',2)
            LegendEntries = {'Observed waveform','Retracking points'};
        end
        plot(1:NrBins,ones(1,NrBins)*CS.Pu(IDXpoi)/NORMfactor(IDXpoi),'k--');
        LegendEntries = horzcat(LegendEntries,{'Pu'});
        axis([0 NrBins 0 max(WD)])
        title('Observed waveform','fontsize',FntSz)
        legend(LegendEntries,'box','off')
        xlabel('Bin','fontsize',FntSz)
        ylabel('Power (Watts)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot coastline
        LOLA = Get_River_Border_data('GSHHS','coastlines','f',DOM(1,:),DOM(2,:));
        subplot(2,3,4),plot(LOLA(:,1),LOLA(:,2),'k-'),hold on
        plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.')
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        axis([minmax(LOLA(:,1)) minmax(LOLA(:,2))])
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)

        %Plot DEM around point of interest including ground track
        if ~isempty(DEM)
        AZIM   = azimuth(CS.GEO.LAT(1:end-1),CS.GEO.LON(1:end-1),CS.GEO.LAT(2:end),CS.GEO.LON(2:end),WGS84);
        MaxDoN = km2deg(20,6371); %Maximum distance-of-nadir [km]
        DEMx   = DEM.GridVectors{1}; DEMy = DEM.GridVectors{2};
        [~,Ix] = min(abs(DEMx-CS.GEO.LON(i))); [~,Iy] = min(abs(DEMy-CS.GEO.LAT(i))); 
        MDidx  = ceil(MaxDoN/mode(diff(DEMx)));
        subplot(2,3,5),imagesc(DEMx(Ix+(-MDidx:MDidx)),DEMy(Iy+(-MDidx:MDidx)),DEM.Values((Iy+(-MDidx:MDidx)),(Ix+(-MDidx:MDidx)))),hold on
        axis square, colorbar, set(gca,'YDir','normal')
        set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String','meters','FontSize',FntSz);
        %Plot coastline
        plot(LOLA(:,1),LOLA(:,2),'-','Color',[.5 .5 .5]),hold on
        %Plot direction of flight
        quiver(CS.GEO.LON(i),CS.GEO.LAT(i),(CS.GEO.LON(i+1)-CS.GEO.LON(i)),(CS.GEO.LAT(i+1)-CS.GEO.LAT(i)),20,'Color','k','LineWidth',2,'MaxHeadSize',0.8);
        %Plot across-track points
        [lat_ac,lon_ac] = reckon(CS.GEO.LAT(i),CS.GEO.LON(i),deg2km(-MaxDoN:mode(diff(DEMx))/2:MaxDoN,6371)*1000,AZIM(i)-90,WGS84);
        plot(lon_ac,lat_ac,'k.-')
        %Plot identified scatterers
        plot(CS.GEO.LON(i),CS.GEO.LAT(i),'o','Color',[.5 .5 .5],'LineWidth',2,'MarkerSize',12)
        title('DEM','fontsize',FntSz)
        xlabel('Longitude (deg)','fontsize',FntSz)
        ylabel('Latitude (deg)','fontsize',FntSz)
        set(gca,'fontsize',FntSz)
        end

        % subplot(2,3,6),plot(CS.GEO.LON(CS.GEO.LON~=0),CS.GEO.LAT(CS.GEO.LON~=0),'.'),hold on
        % plot(CS.GEO.LON(i),CS.GEO.LAT(i),'ro','LineWidth',2)
        % axis([minmax(DEMx(Ix+(-MDidx:MDidx))) minmax(DEMy(Iy+(-MDidx:MDidx)))])
        % plot_google_map
        % xlabel('Longitude (deg)','fontsize',FntSz)
        % ylabel('Latitude (deg)','fontsize',FntSz)
        % set(gca,'fontsize',FntSz)

        %Set background color to white
        set(gcf, 'Color',[1 1 1])
    end
end

%% Save output
DATA            = struct;
IDX             = ~isnan(CS.HEI);
DATA.('TIME')   = CS.TIME(IDX);
DATA.('LAT')    = CS.GEO.LAT(IDX);
DATA.('LON')    = CS.GEO.LON(IDX);
DATA.('HEI')    = single(CS.HEI(IDX));
DATA.('SSHi')   = single(CS.SSHi(IDX));
DATA.('sigma0') = single(CS.sigma0(IDX));
DATA.('SWH')    = single(CS.SWH(IDX));
DATA.('SurfT')  = CS.surf_type(IDX);

end
