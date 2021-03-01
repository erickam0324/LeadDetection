function sigma0 = Compute_Sigma0(CS,DDAcf)

%COMPUTE_SIGMA0 computes the radar backscatter per unit area = sigma naught
%(sigma0). The equations are taken from Dinardo (2016), Guidelines for
%reverting Waveform Power to Sigma Naught for CryoSat-2 in SAR mode. To
%compute the resolution ground-cell in SAR mode (A_SAR), we included the
%term "0.886 .* alpha_earth" as given in the "LOTUS	D2.5 - Snow Depth
%Theoretical Basis Document".

%Input:
%CS:     retracked range to surface reflection point, Transmitted Peak
%        Power, and satellite velocity 

%Output:
%sigma0: sigma naught values in dB

%Define terms to be used in the equations
R            = CS.range;         %Range from Satellite CoM to surface reflection point [m]
Tx_Pwr       = CS.MEA.Tx_Pwr(:); %Transmitted Peak Power [Watt]
L_atm        = 1;                %Two Ways Atmosphere Losses (SO FAR, WE ASSUME THERE ARE NO LOSSES)
L_RX         = 1;                %Receiving Chain (RX) Waveguide Losses (SO FAR, WE ASSUME THERE ARE NO LOSSES)
wf           = 1;                %Footprint widening factor (1 in case of no weighting window application and 1.486*rv in case of Hamming window application on burst data)
R_earth      = 6371000;          %Mean Earth Radius [m]
V_s          = CS.GEO.V.V(:);    %Satellite Along Track Velocity [m/s]

%Compute the resolution ground-cell in SAR mode
alpha_earth  = 1 + (R/R_earth);                              %[Eq. 24]
L_x          = DDAcf.lambda0*R./(2*V_s*DDAcf.tau_B);         %[Eq. 23]
L_y          = sqrt(DDAcf.c*R*DDAcf.PTR_width./alpha_earth); %[Eq. 23]
A_SAR        = (2*L_y).*(wf*L_x);                            %Eq. 22 (Guidelines for reverting Waveform Power to Sigma Naught for CryoSat-2 in SAR mode)
% A_SAR        = (2*L_y).*(wf*L_x) * 0.886 .* alpha_earth;     %Eq. 7 (LOTUS_D25_Snow_Depth_Theoretical_Basis_Document.pdf)

%Determine K defined in Eq. 21
K            = ((4*pi)^3 * R.^4 * L_atm * L_RX) ./ (DDAcf.lambda0^2 * DDAcf.G0^2 * A_SAR);

%Compute sigma naught [Eq. 20]
sigma0       = 10*log10(CS.Pu./Tx_Pwr) + 10*log10(K) + DDAcf.Biases.sigma0;

end
