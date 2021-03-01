classdef RETRACKER
    methods (Static)
        function [fun,x0,lb,ub,XDATA,IDXnr] = Beta5(WD,n,Ypk,Xpk,TrailingEdgeApprox)
            %5-parameter Beta-retracker designed by Martin et al. (Martin,
            %T. V., Zwally, H. J., Brenner, A. C., & Bindschadler, R. A.
            %(1983). Analysis and retracking of continental ice sheet radar
            %altimeter waveforms. Journal of Geophysical Research: Oceans
            %(1978–2012), 88(C3), 1608-1616.)
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -TrailingEdgeApprox = 'linear'/'exponential' trailing edge.
            %
            % Output:
            % -fun     = function handle Beta-5 retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016
            
            %Settings
            nalias  = 12; %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)

            %Initial values
            b1 = mean(WD(1:nalias));             %Thermal noise level
            b2 = Ypk - b1;                       %Return signal amplitude = maximum amplitude over thermal noise
            b4 = Xpk-find(WD > 0.005,1,'first'); %Waveform rise time
            b3 = Xpk - b4/2;                     %Mid-point of leading edge
            b5 = 0;                              %Slope of trailing edge
            x0 = [b1 b2 b3 b4 b5];
    
            %Set input data for objective function
            XDATA = n;
            
            %Set lower bounds, upper bounds, and function handle
            switch TrailingEdgeApprox
                case 'linear'
                    lb  = [0 0.9 find(WD > 0.005,1,'first') 0 -Inf];
                    ub  = [max(WD(1:nalias)) 1 Xpk Xpk 0];
                    fun = @RETRACKER.beta5fun;
                case 'exponential'
                    lb  = [0 0.9 find(WD > 0.005,1,'first') 0 0];
                    ub  = [max(WD(1:nalias)) 1 Xpk Xpk Inf];
                    fun = @RETRACKER.beta5fun_expTE;
                otherwise
                error('TrailingEdgeApprox not reckognized!')
            end
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr = 3;
        end

        function y = beta5fun(x,n)
            %5-parameter Beta-retracker with a linear trailing edge
            %designed by Martin et al. (Martin, T. V., Zwally, H. J.,
            %Brenner, A. C., & Bindschadler, R. A. (1983). Analysis and
            %retracking of continental ice sheet radar altimeter waveforms.
            %Journal of Geophysical Research: Oceans (1978–2012), 88(C3),
            %1608-1616.)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta5;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 5-parameter function
            %
            % Author: D.C. Slobbe, January 2016
            
            B1      = x(1); B2 = x(2); B3 = x(3); B4 = x(4); B5 = x(5);
            Q       = zeros(size(n));
            k       = 0.5;
            IDX     = n >= (B3 + 0.5*B4);
            Q(IDX)  = n(IDX) - (B3 + k*B4);
            y       = B1 + B2.*(1+B5.*Q) .* ( erf((n-B3)/B4/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5;
        end
        
        function y = beta5fun_expTE(x,n)
            %5-parameter Beta-retracker with an exponential trailing edge
            %designed by Deng and Featherstone. (Deng, X., & Featherstone,
            %W. E. (2006). A coastal retracking system for satellite radar
            %altimeter waveforms: Application to ERS‐2 around Australia.
            %Journal of Geophysical Research: Oceans (1978–2012), 111(C6).)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta5;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 5-parameter function
            %
            % Author: D.C. Slobbe, January 2016

            B1      = x(1); B2 = x(2); B3 = x(3); B4 = x(4); B5 = x(5);
            Q       = zeros(size(n));
            k       = 0.5;
            IDX     = n >= (B3 - 2*B4);
            Q(IDX)  = n(IDX) - (B3 + k*B4);
            y       = B1 + B2.*exp(-B5.*Q) .* ( erf((n-B3)/B4/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5;
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = Beta9(WD,n,Ypk,Xpk,Wpk,TrailingEdgeApprox)
            %9-parameter Beta-retracker designed by Martin et al. (Martin,
            %T. V., Zwally, H. J., Brenner, A. C., & Bindschadler, R. A.
            %(1983). Analysis and retracking of continental ice sheet radar
            %altimeter waveforms. Journal of Geophysical Research: Oceans
            %(1978–2012), 88(C3), 1608-1616.)
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -Wpk    = width of peaks;
            % -TrailingEdgeApprox = 'linear'/'exponential' trailing edge.
            %
            % Output:
            % -fun     = function handle Beta-9 retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias  = 12; %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)

            %Initial values peak 1
            b1 = mean(WD(1:nalias));                %Thermal noise level
            b2 = Ypk(1) - b1;                       %Return signal amplitude = maximum amplitude over thermal noise
            b4 = Xpk(1)-find(WD > 0.005,1,'first'); %Waveform rise time
            b3 = Xpk(1) - b4/2 ;                    %Mid-point of leading edge
            b5 = 0;                                 %Slope of trailing edge
            
            %Initial values peak 2
            b6 = Ypk(2) - b1;                       %Return signal amplitude = maximum amplitude over thermal noise
            b8 = Wpk(2);                            %Waveform rise time
            b7 = Xpk(2) - b8/2;                     %Mid-point of leading edge
            b9 = 0;                                 %Slope of trailing edge
            x0 = [b1 b2 b3 b4 b5 b6 b7 b8 b9];
    
            %Set input data for objective function
            XDATA = n;

            %Set lower bounds, upper bounds, and function handle
            switch TrailingEdgeApprox
                case 'linear'
                    lb  = [0 Ypk(1)-(10*b1) find(WD > 0.005,1,'first') 0 -Inf Ypk(2)-(10*b1) Xpk(1) 0 -Inf];
                    ub  = [max(WD(1:nalias)) 1 Xpk(1) Xpk(1) 0 1 Xpk(2) Xpk(2)-Xpk(1) 0];
                    fun = @RETRACKER.beta9fun;
                case 'exponential'
                    lb  = [0 Ypk(1)-(10*b1) find(WD > 0.005,1,'first') 0 0 Ypk(2)-(10*b1) Xpk(1) 0 0];
                    ub  = [max(WD(1:nalias)) 1 Xpk(1) Xpk(1) Inf 1 Xpk(2) Xpk(2)-Xpk(1) Inf];
                    fun = @RETRACKER.beta9fun_expTE;
                otherwise
                error('TrailingEdgeApprox not reckognized!')
            end
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr = 3:4:numel(x0);
        end

        function y = beta9fun(x,n)
            %9-parameter Beta-retracker with a linear trailing edge
            %designed by Martin et al. (Martin, T. V., Zwally, H. J.,
            %Brenner, A. C., & Bindschadler, R. A. (1983). Analysis and
            %retracking of continental ice sheet radar altimeter waveforms.
            %Journal of Geophysical Research: Oceans (1978–2012), 88(C3),
            %1608-1616.)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta9;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 9-parameter function
            %
            % Author: D.C. Slobbe, January 2016
            
            B1      = x(1); B2 = x(2); B3 = x(3); B4 = x(4); B5 = x(5);
            B6      = x(6); B7 = x(7); B8 = x(8); B9 = x(9);
            [Q1,Q2] = deal(zeros(size(n)));
            k       = 0.5;
            IDX     = n >= (B3 + 0.5*B4);
            Q1(IDX) = n(IDX) - (B3 + k*B4);
            IDX     = n >= (B7 + 0.5*B8);
            Q2(IDX) = n(IDX) - (B7 + k*B8);
            y       = B1 + (B2.*(1+B5.*Q1) .* ( erf((n-B3)/B4/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5)...
                         + (B6.*(1+B9.*Q2) .* ( erf((n-B7)/B8/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5);
        end
        
        function y = beta9fun_expTE(x,n)
            %9-parameter Beta-retracker with an exponential trailing edge
            %designed by Deng and Featherstone. (Deng, X., & Featherstone,
            %W. E. (2006). A coastal retracking system for satellite radar
            %altimeter waveforms: Application to ERS‐2 around Australia.
            %Journal of Geophysical Research: Oceans (1978–2012), 111(C6).)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta9;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 9-parameter function
            %
            % Author: D.C. Slobbe, January 2016
            
            B1      = x(1); B2 = x(2); B3 = x(3); B4 = x(4); B5 = x(5);
            B6      = x(6); B7 = x(7); B8 = x(8); B9 = x(9);
            [Q1,Q2] = deal(zeros(size(n)));
            k       = 0.5;
            IDX     = n >= (B3 - 2*B4);
            Q1(IDX) = n(IDX) - (B3 + k*B4);
            IDX     = n >= (B7 - 2*B8);
            Q2(IDX) = n(IDX) - (B7 + k*B8);
            y       = B1 + (B2.*exp(-B5.*Q1) .* ( erf((n-B3)/B4/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5)...
                         + (B6.*exp(-B9.*Q2) .* ( erf((n-B7)/B8/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5);
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = BetaX(WD,n,Ypk,Xpk,Wpk,TrailingEdgeApprox)
            %X-parameter Beta-retracker designed by Martin et al. (Martin,
            %T. V., Zwally, H. J., Brenner, A. C., & Bindschadler, R. A.
            %(1983). Analysis and retracking of continental ice sheet radar
            %altimeter waveforms. Journal of Geophysical Research: Oceans
            %(1978–2012), 88(C3), 1608-1616.)
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -Wpk    = width of peaks;
            % -TrailingEdgeApprox = 'linear'/'exponential' trailing edge.
            %
            % Output:
            % -fun     = function handle Beta-X retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias  = 12; %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)

            %Initial values
            NumP   = numel(Ypk);
            b1     = mean(WD(1:nalias));                %Thermal noise level
            %[Return signal amplitude = maximum amplitude over thermal noise; Mid-point of leading edge; Waveform rise time; Slope of trailing edge]
            b      = [Ypk' - b1;Xpk' - 0.5*Wpk';Wpk';zeros(1,NumP)];
            b(3,1) = Xpk(1)-find(WD > 10*b1,1,'first');
            b(3,1) = max([b(3,1);0]);
            b(2,:) = Xpk' - b(3,:)./2;
            x0     = [b1 b(:)'];

            %Set input data for objective function
            XDATA = n;

            %Set lower bounds, upper bounds, and function handle
            switch TrailingEdgeApprox
                case 'linear'
                    lb  = [Ypk' - (10*b1);[find(WD > 0.005,1,'first') Xpk(1:end-1)'];zeros(1,NumP);-Inf(1,NumP)];
                    lb  = [0 lb(:)'];
                    ub  = [ones(1,NumP);Xpk';[Xpk(1) Xpk(2:end)'-Xpk(1:end-1)'];zeros(1,NumP)];
                    ub  = [max(WD(1:nalias)) ub(:)'];
                    fun = @RETRACKER.betaXfun;
                case 'exponential'
                    lb  = [Ypk' - (10*b1);[find(WD > 0.005,1,'first') Xpk(1:end-1)'];zeros(2,NumP)];
                    lb  = [0 lb(:)'];
                    ub  = [ones(1,NumP);Xpk';[Xpk(1) Xpk(2:end)'-Xpk(1:end-1)'];Inf(1,NumP)];
                    ub  = [max(WD(1:nalias)) ub(:)'];
                    fun = @RETRACKER.betaXfun_expTE;
                otherwise
                error('TrailingEdgeApprox not reckognized!')
            end
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr = 3:4:numel(x0);
        end

        function y = betaXfun(x,n)
            %9-parameter Beta-retracker with a linear trailing edge
            %designed by Martin et al. (Martin, T. V., Zwally, H. J.,
            %Brenner, A. C., & Bindschadler, R. A. (1983). Analysis and
            %retracking of continental ice sheet radar altimeter waveforms.
            %Journal of Geophysical Research: Oceans (1978–2012), 88(C3),
            %1608-1616.)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta9;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 9-parameter function
            %
            % Author: D.C. Slobbe, January 2016
            
            k = 0.5;
            y = ones(size(n))*x(1);
            for i=1:(numel(x)-1)/4
                Q      = zeros(size(n));
                IDX    = n >= (x(3+((i-1)*4)) + 0.5*x(4+((i-1)*4)));
                Q(IDX) = n(IDX) - (x(3+((i-1)*4)) + k*x(4+((i-1)*4)));
                y      = y + (x(2+((i-1)*4)).*(1+x(5+((i-1)*4)).*Q) .* ( erf((n-x(3+((i-1)*4)))/x(4+((i-1)*4))/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5);
            end
        end
        
        function y = betaXfun_expTE(x,n)
            %9-parameter Beta-retracker with an exponential trailing edge
            %designed by Deng and Featherstone. (Deng, X., & Featherstone,
            %W. E. (2006). A coastal retracking system for satellite radar
            %altimeter waveforms: Application to ERS‐2 around Australia.
            %Journal of Geophysical Research: Oceans (1978–2012), 111(C6).)
            %
            %Input: 
            % -x      = vector with Beta1 to Beta9;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the 9-parameter function
            %
            % Author: D.C. Slobbe, January 2016
            
            k = 0.5;
            y = ones(size(n))*x(1);
            for i=1:(numel(x)-1)/4
                Q      = zeros(size(n));
                IDX    = n >= (x(3+((i-1)*4)) - 2*x(4+((i-1)*4)));
                Q(IDX) = n(IDX) - (x(3+((i-1)*4)) + k*x(4+((i-1)*4)));
                y      = y + (x(2+((i-1)*4)).*exp(-x(5+((i-1)*4)).*Q) .* ( erf((n-x(3+((i-1)*4)))/x(4+((i-1)*4))/sqrt(2)) - erf(-inf/sqrt(2)) )/2^1.5);
            end
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = Brown(WD,n,Ypk,Xpk,h,dt)
            %Brown's theoretical ocean model after the equations given by
            %Passaro et al. (Passaro, M., Cipollini, P., Vignudelli, S.,
            %Quartly, G. D., & Snaith, H. M. (2014). ALES: A multi-mission
            %adaptive subwaveform retracker for coastal and open ocean
            %altimetry. Remote Sensing of Environment, 145, 173-189.). To
            %understand his derivation, use Amarouche, L., Thibaut, P.,
            %Zanife, O. Z., Dumont, J. P., Vincent, P., & Steunou, N.
            %(2004). Improving the Jason-1 ground retracking to better
            %account for attitude effects. Marine Geodesy, 27(1-2),
            %171-197.
            %
            %Originally developed by Brown, G. S. (1977). The average
            %impulse response of a rough surface and its applications.
            %Antennas and Propagation, IEEE Transactions on, 25(1), 67-74.
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -h      = satellite altitude [m];
            % -dt     = sampling interval / compressed pulse length [ns]
            %
            % Output:
            % -fun     = function handle Brown's retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias  = 12;        %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)
            c       = 299792458; %Light velocity [m/s]

            %Initial values
            Tn       = mean(WD(1:nalias));               %Thermal noise level
            Pu       = Ypk(1) - Tn;                      %Amplitude of the signal (relates to the backscatter coefficient sigma0)
            tau      = RETRACKER.Threshold(WD,n',WD)*dt; %Epoch or time delay i.e. the position of the waveform in the analysis window, w.r.t. nominal tracking reference point [ns]
            SWH      = 2;
            sigma_s  = SWH/(2*c);                        %Leading edge slope (relates to the significant wave height)
            xi       = 0;                                %Off-nadir mispointing angle
            x0       = [Tn Pu tau sigma_s xi];
            
            %Set input data for objective function
            XDATA    = [n;h;dt/1E9];

            %Set lower bounds, upper bounds, and function handle
            lb  = [0 1-Tn find(WD > 0.005,1,'first')*dt 0 0];
            ub  = [max(WD(1:nalias)) 1 Xpk*dt 50/(2*c) Inf];
            fun = @RETRACKER.Brownfun;
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr = 3;
        end
        
        function y = Brownfun(x,n)
            %Brown's theoretical ocean model after the equations given by
            %Passaro et al. (Passaro, M., Cipollini, P., Vignudelli, S.,
            %Quartly, G. D., & Snaith, H. M. (2014). ALES: A multi-mission
            %adaptive subwaveform retracker for coastal and open ocean
            %altimetry. Remote Sensing of Environment, 145, 173-189.)   
            %
            %Input: 
            % -x      = vector with unknowns;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the Brown's model
            %
            % Author: D.C. Slobbe, January 2016
            
            %Copy satellite altitude to variable "h" and sampling interval
            %/ compressed pulse length to variable "rt"
            h        = n(end-1);
            rt       = n(end);
            n        = n(1:end-2);
            t        = n*rt;
            
            c        = 299792458;       %Light velocity [m/s]
            Re       = 6371000;         %Mean Earth radius [m]
            sigma_p  = 0.513*rt;        %Width of the radar point target response function (Hayne, 1980)
            theta0   = deg2rad(1.06);   %Antenna beam width (along-track) [radians]
%             theta0  = deg2rad(1.1992); %Antenna beam width (across-track) [radians]

            Tn       = x(1);     %Thermal noise level
            Pu       = x(2);     %Amplitude of the signal (relates to the backscatter coefficient sigma0)
            tau      = x(3)/1E9; %Epoch or time delay i.e. the position of the waveform in the analysis window, w.r.t. nominal tracking reference point
            sigma_s  = x(4);     %Leading edge slope (relates to the significant wave height)
            xi       = x(5);     %Off-nadir mispointing angle
            
            gamma    = sin(theta0)^2 / (2*log(2));
            a_xi     = exp((-4*sin(xi)^2)/gamma);

            b_xi     = cos(2*xi) - sin(2*xi)^2/gamma;
            a        = (4*c) / (gamma*h*(1 + h/Re));
            c_xi     = b_xi*a;
            sigma_c  = sqrt(sigma_p^2 + sigma_s^2);
            v        = c_xi .* (t-tau-(c_xi*sigma_c^2)/2);
            
            u        = (t-tau-c_xi*sigma_c^2) / (sqrt(2)*sigma_c);
            
            y        = Tn + (a_xi .* Pu  .* ((1 + erf(u))/2) .* exp(-v));
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = BrownHayne(WD,n,Ypk,Xpk,h,dt)
            %Brown-Hayne Theoretical Ocean Model after the equations given
            %by Gommenginger et al. (Gommenginger, C., Thibaut, P.,
            %Fenoglio-Marc, L., Quartly, G., Deng, X., Gómez-Enri, J., ...
            %& Gao, Y. (2011). Retracking altimeter waveforms near the
            %coasts. In Coastal altimetry (pp. 61-101). Springer Berlin
            %Heidelberg.)
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -h      = satellite altitude [m];
            % -dt     = sampling interval / compressed pulse length [ns].
            %
            % Output:
            % -fun     = function handle Brown-Hayne retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias  = 12;        %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)
            c       = 299792458; %Light velocity [m/s]

            %Initial values
            Tn       = mean(WD(1:nalias));               %Thermal noise level
            Pu       = Ypk(1) - Tn;                      %Amplitude of the signal (relates to the backscatter coefficient sigma0)
            tau      = RETRACKER.Threshold(WD,n',WD)*dt; %Epoch or time delay i.e. the position of the waveform in the analysis window, w.r.t. nominal tracking reference point [ns]
            SWH      = 2;
            sigma_s  = SWH/(2*c);                        %Leading edge slope (relates to the significant wave height)
            xi       = 0;                                %Off-nadir mispointing angle
            lambda_s = 0;                                %Skewness
            x0       = [Tn Pu tau sigma_s xi lambda_s];

            %Set input data for objective function
            XDATA    = [n;h;dt/1E9];

            %Set lower bounds, upper bounds, and function handle
            lb  = [0 1-Tn find(WD > 0.005,1,'first')*dt 0 0 -Inf];
            ub  = [max(WD(1:nalias)) 1 Xpk*dt 50/(2*c) Inf Inf];
            fun = @RETRACKER.BrownHaynefun;
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr = 3;
        end
        
        function y = BrownHaynefun(x,n)
            %Brown-Hayne Theoretical Ocean Model after the equations given
            %by Gommenginger et al. (Gommenginger, C., Thibaut, P.,
            %Fenoglio-Marc, L., Quartly, G., Deng, X., Gómez-Enri, J., ...
            %& Gao, Y. (2011). Retracking altimeter waveforms near the
            %coasts. In Coastal altimetry (pp. 61-101). Springer Berlin
            %Heidelberg.)
            %
            %Input: 
            % -x      = vector with unknowns;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the Brown-Hayne model
            %
            % Author: D.C. Slobbe, January 2016
            
            %Copy satellite altitude to variable "h" and sampling interval
            %/ compressed pulse length to variable "rt"
            h        = n(end-1);
            rt       = n(end);
            n        = n(1:end-2);
            t        = n*rt;

            c        = 299792458;       %Light velocity [m/s]
            Re       = 6371000;         %Mean Earth radius [m]
            sigma_p  = 0.513*rt;        %Width of the radar point target response function (Hayne, 1980)
            theta0   = deg2rad(1.06);   %Antenna beam width (along-track) [radians]
%             theta0  = deg2rad(1.1992); %Antenna beam width (across-track) [radians]

            Tn       = x(1);     %Thermal noise level
            Pu       = x(2);     %Amplitude of the signal (relates to the backscatter coefficient sigma0)
            tau      = x(3)/1E9; %Epoch or time delay i.e. the position of the waveform in the analysis window, w.r.t. nominal tracking reference point
            sigma_s  = x(4);     %Leading edge slope (relates to the significant wave height)
            xi       = x(5);     %Off-nadir mispointing angle
            lambda_s = x(6);     %Skewness
            
            gamma    = sin(theta0)^2 / (2*log(2));
            a_xi     = exp((-4*sin(xi)^2)/gamma);

            b_xi     = cos(2*xi) - sin(2*xi)^2/gamma;
            a        = (4*c) / (gamma*h*(1 + h/Re));
            c_xi     = b_xi*a;
            sigma_c  = sqrt(sigma_p^2 + sigma_s^2);
            v        = c_xi .* (t-tau-(c_xi*sigma_c^2)/2);
            
            u        = (t-tau-c_xi*sigma_c^2) / (sqrt(2)*sigma_c);
            
            y        = Tn + a_xi .* (Pu/2) .* exp(-v) .* ( (1 + erf(u)) + ( ((lambda_s/6)*(sigma_s/sigma_c)^2) .* ...
                ( (1 + erf(u))*c_xi^3*sigma_c^3 - ((sqrt(2)/sqrt(pi)) .* (2*u.^2 + 3*sqrt(2)*c_xi*sigma_c.*u + ...
                3*c_xi^2*sigma_c^2 - 1) .* exp(-u.^2) ) ) ) );
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = D2P(WD,n,Ypk,Xpk,Wpk)
            %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
            %designed by Giles et al. (Giles, K. A., Laxon, S. W., Wingham,
            %D. J., Wallis, D. W., Krabill, W. B., Leuschen, C. J., ... &
            %Raney, R. K. (2007). Combined airborne laser and radar
            %altimeter measurements over the Fram Strait in May 2002.
            %Remote Sensing of Environment, 111(2), 182-194.
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -Wpk    = width of peaks.
            %
            % Output:
            % -fun     = function handle Brown-Hayne retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias = 12; %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)

            %Initial values
            NumP   = numel(Ypk);
            b1     = mean(WD(1:nalias)); %Thermal noise level
            k      = 0.1;                %decay of trailing edge
            %[max amplitude of echo; time of max amplitude; dist(end of leading edge) - dist(start of trailing edge); decay of trailing edge; width of leading edge]
            x0     = reshape([Ypk Xpk 0.5*Wpk ones(NumP,1)*k Wpk]',NumP*5,1);
            
            %Set input data for objective function
            XDATA  = n;

            %Set lower bounds, upper bounds, and function handle
            if NumP == 1
                lb = [Ypk' - (10*b1);find(WD > 0.005,1,'first');zeros(3,NumP)];
                ub = [1;Xpk';min([Xpk,max(n)-Xpk]);Inf;min([Xpk,max(n)-Xpk])];
            else
                lb = [Ypk' - (10*b1);[find(WD > 0.005,1,'first') Xpk(1:end-1)'];zeros(3,NumP)];
                ub = [ones(1,NumP);Xpk';0.5*[diff(Xpk'),diff(Xpk(end-1:end))];Inf(1,NumP);0.5*[diff(Xpk'),diff(Xpk(end-1:end))]];
            end
            fun    = @RETRACKER.D2Pfun;
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr  = 2:5:numel(x0);
        end
        
        function y = D2Pfun(x,n)
            %D2P (Delay/Doppler Phase-monopulse Radar Altimeter) retracker
            %designed by Giles et al. (Giles, K. A., Laxon, S. W., Wingham,
            %D. J., Wallis, D. W., Krabill, W. B., Leuschen, C. J., ... &
            %Raney, R. K. (2007). Combined airborne laser and radar
            %altimeter measurements over the Fram Strait in May 2002.
            %Remote Sensing of Environment, 111(2), 182-194.
            %
            %Input: 
            % -x      = vector with unknowns;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the D2P function
            %
            % Author: D.C. Slobbe, January 2016

            npeaks = numel(x)/5;
            x      = reshape(x,5,npeaks);
            [f,y]  = deal(zeros(size(n)));
            for i = 1:npeaks
                a       = x(1,i); % maximal return amplitude (in this case, amplitude of peak)
                t0      = x(2,i); % time of maximal return amplitude
                tb      = x(3,i); % dist(end of leading edge) - dist(start of trailing edge)
                k       = x(4,i); % decay of trailing edge
                sigma   = x(5,i); % width of leading edge
                a3      = ((-sqrt(k*tb)*2*sigma/tb) + 2) / ((tb^2 + 2*t0^2 + 2*t0*tb) * (3 - 2*sigma));
                a2      = (sqrt(k*tb)/tb^2) - (1/(sigma*tb)) - (a3*tb);
                IDX1    = n < t0;
                f(IDX1) = (n(IDX1) - t0)./sigma;
                IDX2    = n >= t0 & n < (t0+tb);
                f(IDX2) = (a3.*(n(IDX2) - t0).^3) + (a2.*(n(IDX2) - t0).^2) + ((n(IDX2) - t0)/sigma);
                IDX3    = n >= (t0+tb);
                f(IDX3) = sqrt(k.*(n(IDX3) - t0));
                y       = y + a.*exp(-f.^2);
            end
        end
        
        function [fun,x0,lb,ub,XDATA,IDXnr] = FunctionFit(WD,n,Ypk,Xpk,Wpk)
            %"Function Fit" retracking is ALT_RET_ICE_03. (Surface
            %Topography Mission (STM) SRAL/MWR L2 Algorithms Definition,
            %Accuracy and Specification [SD-03] [SD-07])
            %
            %Input: 
            % -WD     = full waveform or sub-waveform;
            % -n      = bin indices of (sub-)waveform;
            % -Ypk    = waveform amplitude of peak that needs to be retracked;
            % -Xpk    = bin index of peak that needs to be retracked;
            % -Wpk    = width of peaks.
            %
            % Output:
            % -fun     = function handle Brown-Hayne retracker;
            % -x0      = vector with initial values of unknown parameters;
            % -lb      = vector of lower bounds;
            % -ub      = vector of upper bounds;
            % -XDATA   = input data for objective function.
            % -IDXnr   = indices in x refering to retracking location
            %
            % Author: D.C. Slobbe, February 2016

            %Settings
            nalias = 12; %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)

            %Initial values
            NumP   = numel(Ypk);
            b1     = mean(WD(1:nalias)); %Thermal noise level
            k      = 0.1;                %decay of trailing edge
            %[max amplitude of echo; time of max amplitude; decay of trailing edge; width of leading edge]
            x0     = reshape([Ypk Xpk ones(NumP,1)*k Wpk]',NumP*4,1);
            
            %Set input data for objective function
            XDATA  = n;

            %Set lower bounds, upper bounds, and function handle
            if NumP == 1
                lb = [Ypk - (10*b1);find(WD > 0.005,1,'first');zeros(2,NumP)];
                ub = [1;Xpk;1;min([Xpk,max(n)-Xpk])];
            else
                lb = [Ypk' - (10*b1);[find(WD > 0.005,1,'first') Xpk(1:end-1)'];zeros(2,NumP)];
                ub = [ones(1,NumP);Xpk';ones(1,NumP);0.5*[diff(Xpk'),diff(Xpk(end-1:end))]];
            end
            fun    = @RETRACKER.FunctionFitfun;
            
            %Set indices of vector with estimated unknowns refering to retracking location
            IDXnr  = 2:4:numel(x0);
        end
        
        function y = FunctionFitfun(x,n)
            %"Function Fit" retracking is ALT_RET_ICE_03. (Surface
            %Topography Mission (STM) SRAL/MWR L2 Algorithms Definition,
            %Accuracy and Specification [SD-03] [SD-07])
            %
            %Input: 
            % -x      = vector with unknowns;
            % -n      = bin indices of (sub-)waveform.
            %
            % Output:
            % -y      = the "Function Fit" retracking function
            %
            % Author: D.C. Slobbe, January 2016

            npeaks = numel(x)/4;
            x      = reshape(x,4,npeaks);
            [f,y]  = deal(zeros(size(n)));
            for i = 1:npeaks
                a       = x(1,i); % maximal return amplitude (in this case, amplitude of peak)
                t0      = x(2,i); % time of maximal return amplitude
                k       = x(3,i); % decay of trailing edge
                sigma   = x(4,i); % width of leading edge
                tb      = k*sigma^2;
                c       = sqrt(k*tb);
                a3      = -(3*k*sigma - 2*c)/(2*sigma*tb^2*c);
                a2      = (-5*k*sigma + 4*c)/(2*sigma*tb*c);
                IDX1    = n < t0;
                f(IDX1) = (n(IDX1) - t0)./sigma;
                IDX2    = n >= t0 & n < (t0+tb);
                f(IDX2) = (a3.*(n(IDX2) - t0).^3) - (a2.*(n(IDX2) - t0).^2) + ((n(IDX2) - t0)/sigma);
                IDX3    = n >= (t0+tb);
                f(IDX3) = sqrt(k.*(n(IDX3) - t0));
                y       = y + a.*exp(-f.^2);
            end
        end
        
        function [n_ret,A] = OCOG(WD,n)
            %Offset Centre Of Gravity retracker designed by Wingham et al. 
            %(Wingham, D. J., Rapley, C. G., & Griffiths, H. (1986). New
            %techniques in satellite altimeter tracking systems. In ESA
            %Proceedings of the 1986 International Geoscience and Remote
            %Sensing Symposium(IGARSS'86) on Remote Sensing: Today's
            %Solutions for Tomorrow's Information Needs, (Vol. 3).)
            %
            %Input: 
            % -n      = bin indices of (sub-)waveform;
            % -WD     = full waveform or sub-waveform.
            %
            % Output:
            % -n_ret  = retracking location [bins];
            % -A      = amplitude
            %
            % Author: D.C. Slobbe, January 2016

            %Amplitude (see Eq. 4.1, Vignudelli et al. Coastal Altimetry)
            A     = sqrt(sum(WD.^4)./sum(WD.^2));

            %Pulse width (see Eq. 4.2, Vignudelli et al. Coastal Altimetry)
            W     = sum(WD.^2).^2./sum(WD.^4);
            
            %Centre Of Gravity (see Eq. 4.3, Vignudelli et al. Coastal Altimetry)
            COG   = diag(n*(WD.^2))'./sum(WD.^2);
            
            %Leading Edge Position (see Eq. 4.4, Vignudelli et al. Coastal Altimetry)
            n_ret = COG-(W/2);
        end
        
        function [x0,lb,ub,XDATA,IDXnr,BeamIDX] = SAMOSA(WD,n,DDAcf,TN,Pu0,t0_0,SWH0,nu0,Re,h,vt,OrbSlp,theta_p,theta_r,SmodID,EstPar,BAstck)
            %3/4-parameter SAMOSA2/SAMOSA3 retrackers described by Ray et
            %al. (2015) & Dinardo et al. (2018).
            %
            %Ray, C., Martin-Puig, C., Clarizia, M. P., Ruffini, G.,
            %Dinardo, S., Gommenginger, C., & Benveniste, J. (2015). SAR
            %altimeter backscattered waveform model. IEEE Transactions on
            %Geoscience and Remote Sensing, 53(2), 911-919.
            %
            %Dinaro, S., Fenoglio-Marc, L., Buchhaupt, C. Becker, M.,
            %Scharroo, R., Fernandes, M.J., & Benveniste, J. (2018).
            %Coastal and PLRM altimetry in German Bight and West Baltic
            %Sea. Advances in Space Research, 62, 1371-1404.
            %
            %Input:
            % -WD      = full waveform
            % -n       = bin indices of waveform
            % -TN      = thermal noise floor
            % -Pu0     = initial value Pu
            % -t0_0    = initial value epoch (t0) with respect to tracker [ns] 
            % -SWH0    = initial value SWH [m]
            % -nu0     = initial value inverse of the mean-square slope of the sea surface
            % -Re      = local radius of curvature of the Earth's surface
            % -h       = presumed distance between instrument and sea surface as estimated by the altimeter's onboard tracker [m]
            % -vt      = tangential velocity [m/s]
            % -OrbSlp  = slope of orbit
            % -theta_p = pitch angle [rad]
            % -theta_r = roll angle [rad]
            % -SmodID  = [1/0] to apply, respectively, SAMOS2 or SAMOSA3
            % -EstPar  = [3/4] to estimate, respectively, [Pu,t0,SWH] or [Pu,t0,nu0]. If other than 3/4, all four parameters will be estimated
            % -BAstck  = Stack with Doppler beam angles
            %
            %Output:
            % -x0      = vector with initial values
            % -lb      = vector with lower bound values
            % -ub      = vector with upper bound values
            % -XDATA   = bin indices WF, TN, SWH, loc. radius of curvature 
            %            of the Earth's surface, dist. satellite-surface,
            %            tangential velocity, slope of orbit, pitch & roll
            %            angles, SmodID, and EstPar
            % -IDXnr   = index of parameter t in vector with unknowns
            % 
            % Author: Marcel Kleinherenbrink (original version, September
            % 2018). D.C. Slobbe (December 2018)
            
            %Settings
            Ns     = length(WD);  %Nr of bins/samples in any waveform
            os_ZP  = Ns/DDAcf.Np; %Zero-Padding Oversampling Factor
            RefBin = (DDAcf.RefBin-1)*os_ZP + 1; %Reference bin
            t_res  = 1E9/DDAcf.B; %Time resolution [ns]
            t_s    = t_res/os_ZP; %Waveform sampling interval [ns] (~= t_res because of zero-padding)
            alpha  = 1+h/Re;      %Geometric parameter
            
            %Initial values
            if isempty(Pu0) || isnan(Pu0),   Pu0 = 1; end
            if isempty(t0_0)
                [~,I] = max(WD);
                t0_0  = -(RefBin-I)*t_s;
            end
            if isempty(SWH0) || isnan(SWH0), SWH0 = 2; end
            if isempty(nu0) || isnan(nu0),   nu0 = 0; end
            x0 = [Pu0,t0_0,SWH0,nu0];

            %ML2. Calculate idealised Doppler beam angles if not available from the SAR L1B input file
            if isempty(BAstck)
                l       = -106:106;
                dTheta  = (vt*DDAcf.BRI) / (h*alpha);        %Eq. 85
                Theta1  = pi/2 + dTheta * l(1);              %Eq. 86
                Theta2  = pi/2 + dTheta * l(end);            %Eq. 87
                BAstck  = Theta1:dTheta:Theta2;              %Eq. 88
            end
            
            %ML3. Calculate the vector of Doppler frequency of the Doppler beams used for multi-looking
            DFstck  = (2*vt/DDAcf.lambda0) *cos(BAstck);     %Eq. 90
 
            %ML4. Sub-setting the Stack
            dfa     = DDAcf.fp/DDAcf.Nb;                     %Eq. 61
            BeamIDX = double(unique(int32(DFstck/dfa)));
            
            %Set input data for objective function
            XDATA = [n;TN;NaN;Re;h;vt;OrbSlp;theta_p;theta_r;SmodID;EstPar];

            %Set lower and upper bounds
            lb    = [0 -RefBin*t_s -1 0];
            [~,I] = max(WD);
            ub    = [40 (I-RefBin)*t_s 20 Inf];
            
            %Set index of parameter refering to the retracking location in vector with unknowns 
            IDXnr = 2;
        end
        
        function y = SAMOSAfun(x,n,DDAcf,LUT_B14,LUT_Bm14,LUT_B34,LUT_Bm34,l)
            %3/4-parameter SAMOSA2/SAMOSA3 retrackers described by Ray et
            %al. (2015) & Dinardo et al. (2018).
            %
            %Ray, C., Martin-Puig, C., Clarizia, M. P., Ruffini, G.,
            %Dinardo, S., Gommenginger, C., & Benveniste, J. (2015). SAR
            %altimeter backscattered waveform model. IEEE Transactions on
            %Geoscience and Remote Sensing, 53(2), 911-919.
            %
            %Dinaro, S., Fenoglio-Marc, L., Buchhaupt, C. Becker, M.,
            %Scharroo, R., Fernandes, M.J., & Benveniste, J. (2018).
            %Coastal and PLRM altimetry in German Bight and West Baltic
            %Sea. Advances in Space Research, 62, 1371-1404.
            %
            %Input: 
            % -x        = vector with unknowns;
            % -n        = bin indices WF, TN, SWH, loc. radius of curvature 
            %             of the Earth's surface, dist. satellite-surface,
            %             tangential velocity, slope of orbit, pitch & roll
            %             angles, SmodID, and EstPar
            % -DDAcf    = DDA configuration file
            % -LUT_B14  = look-up-table besseli(1/4,LUT_x,1)
            % -LUT_Bm14 = look-up-table besseli(-1/4,LUT_x,1)
            % -LUT_B34  = look-up-table besseli(3/4,LUT_x,1)
            % -LUT_Bm34 = look-up-table besseli(-3/4,LUT_x,1)
            % -l        = vector with Doppler beam indices
            %
            % Output:
            % -y        = the SAMOSA(2/3)(FF) model
            %
            % Author: Marcel Kleinherenbrink (original version, September
            % 2018). D.C. Slobbe (December 2018)
            
            %Unpack vector n
            EstPar    = n(end);
            SmodID    = n(end-1);
            theta_r   = n(end-2);
            theta_p   = n(end-3);
            Orb_slope = n(end-4);
            vt        = n(end-5);
            h         = n(end-6);
            Re        = n(end-7);
            DUMx      = n(end-8);
            TN        = n(end-9);
            n         = n(1:end-10);
            
            %Unpack vector with unknowns. Convert X2 to seconds
            if EstPar == 3
                X1 = x(1); X2=x(2)/1E9; X3=x(3); X4=0;
            elseif EstPar == 4
                X1 = x(1); X2=x(2)/1E9; X3=DUMx; X4=x(end)*1E5;
            else
                X1 = x(1); X2=x(2)/1E9; X3=x(3); X4=x(4)*1E5;
            end

            %Preliminaries
            Ns        = length(n);                                     %Nr of bins/samples in any waveform (after zero-padding)
            os_ZP     = Ns/DDAcf.Np;                                   %Zero-Padding Oversampling Factor
            RefBin    = (DDAcf.RefBin-1)*os_ZP + 1;                    %Reference bin
            alpha     = 1+h/Re;                                        %Geometric parameter
            t_res     = 1/DDAcf.B;                                     %Time resolution [s]
            t_s       = t_res/os_ZP;                                   %Waveform sampling interval [s] (~= t_res because of zero-padding)
            Lx        = (DDAcf.c*h*DDAcf.fp)/(2*vt*DDAcf.fc*DDAcf.Nb); %Along-track resolution (Ray et al. (2015), Eq. 14)
            Ly        = sqrt((DDAcf.c*h)/(alpha*DDAcf.B));             %Ray et al. (2015), Eq. 21
            Lz        = DDAcf.c/(2*DDAcf.B);                           %Range resolution (Ray et al. (2015), Eq. 21)
            gamma_x   = 8*log(2)/DDAcf.theta_x^2;                      %Along-track antenna gain parameter (Dinardo et al. (2018), Eq. 27)
            gamma_y   = 8*log(2)/DDAcf.theta_y^2;                      %Cross-track antenna gain parameter (Dinardo et al. (2018), Eq. 27)
            L_gamma   = alpha*h/(2*gamma_y);                           %(Dinardo et al. (2018), Eq. 27)
            %To mitigate the effect of the model's approximation of the
            %squared PTR with a Gaussian curve (Fenoglio-Marc et al.,
            %2015a), a dynamic alpha_p value is used depending on the SWH
            %for waveforms not in the coastal zone. Inside, the coastal
            %zone, Dinardo et al. (2018) use a constant value of 0.55;
%             if EstPar == 4
%                 alpha_p = 0.552;
%             else
                %The coefficients of the polynomial p of degree 21 are
                %obtained using the polyfit routine with as input the
                %look-up-table alphaPTR_table_05_Aug_2014.txt provided by
                %Dinardo
                % LUT      = load('alphaPTR_table_05_Aug_2014.txt');
                % [p,S,mu] = polyfit(LUT(:,1),LUT(:,2),21);
                p  = [-3.32997745E-05 4.22960667E-05 4.50359245E-04 -5.36103833E-04 ...
                    -2.62882257E-03 2.90968939E-03 8.61154298E-03 -8.75658605E-03 ...
                    -1.73230951E-02 1.59306322E-02 2.19834910E-02 -1.78660221E-02 ...
                    -1.75148616E-02 1.22405603E-02 8.21268631E-03 -4.50013043E-03 ...
                    -2.62173345E-03 1.77289338E-03 -1.35818266E-03 2.00676416E-03 ...
                    1.90840697E-01 7.18194984E-01];
                mu = [1.00000000E+01 5.77783264E+00];
                if X3 > 20
                    alpha_p = polyval(p,20,[],mu);
                elseif X3 < 0
                    alpha_p = polyval(p,0,[],mu);
                else
                    alpha_p = polyval(p,X3,[],mu);
                end
%             end
            sigma_p   = alpha_p*t_res;               %Standard deviation of the Gaussian function that models the PTR, related to the effective pulse length by sigma_p = alpha_p * tau_p (Amarouche et al. 2004). tau_p = 3.125 ns
            theta_lim = Lz/(alpha*Lx);               %Max along-track look angle variation within a strip (Dinardo et al. (2018), Eq. 18)
            sigma_z   = X3/4;                        %(Dinardo et al. (2018), pp. 1381)
            sigma_s   = sigma_z/Lz;                  %(Ray et al. (2015), Eq. 37)
            
            %Get SAMOSA analytical model
            t = n*t_s-t_s*RefBin;
            
            %Apply slope effect compensation in the echo model (Gommenginger et al. (2017), Eqs. 33,75)
            l       = l-(Orb_slope * h / (alpha * Lx));
            
            %Look angle (Dinardo et al. (2018), Eq. 16)
            theta_l = Lx/h*l;
            
            %gamma_0 (Dinardo et al. (2018), Eq. 25)
            gamma_0 = exp(-gamma_y*theta_r^2-X4*theta_l.^2-gamma_x*(theta_l-theta_p).^2)...
                .*exp(-(gamma_y+X4)*(DDAcf.c*(t-X2)/(alpha*h))).*cosh(2*gamma_y*theta_r*(DDAcf.c*(t-X2)/(alpha*h)).^0.5);
            gamma_0(t<=X2,:)=repmat(exp(-gamma_y*theta_r^2-X4*theta_l.^2-gamma_x*(theta_l-theta_p).^2),sum(t<=X2),1);
            
            %T_k (Dinardo et al. (2018), Eq. 26)
            T_k        = (1+X4/gamma_y)-theta_r./(DDAcf.c*(t-X2)/(alpha*h)).^0.5.*tanh(2*gamma_y*theta_r*(DDAcf.c*(t-X2)/(alpha*h)).^0.5);
            T_k(t<=X2) = (1+X4/gamma_y)-2*gamma_y*theta_r^2;
            
            %sigma_c (Dinardo et al. (2018), Eq. 20)
            % sigma_c = (sigma_p^2*(1+(theta_l/theta_lim).^2)+X3/abs(X3)*(abs(X3)/(2*DDAcf.c))^2).^0.5;
            sigma_c = sqrt(sigma_p^2*(1+(Lx/Lz)^2*(alpha*theta_l).^2)+(2*sigma_z/DDAcf.c)^2);
            
            %g_l (Dinardo et al. (2018), Eqs. 17 and 20)
            % g_l = 1./sqrt(alpha_p^2+4*alpha_p^2*(Lx/Ly)^4*l.^2+(sigma_z/Lz)^2);
            % g_l = 1./(DDAcf.B*sqrt(sigma_p^2+4*sigma_p^2*(Lx/Ly)^4*l.^2+(2*sigma_z/DDAcf.c)^2));
            g_l = 1./(DDAcf.B*sqrt(sigma_p^2*(1+(Lx/Lz)^2*(alpha*theta_l).^2)+(2*sigma_z/DDAcf.c)^2));
            % g_l = 1./(sigma_c*DDAcf.B);

            %Use LUTs to evaluate the modified bessel functions of the
            %first kind
            xx                  = ((t-X2)./(2*sigma_c)).^2;
            [B14,Bm14,B34,Bm34] = deal(nan(size(xx)));
            B14(:)              = LUT_B14(xx(:));
            Bm14(:)             = LUT_Bm14(xx(:));
            B34(:)              = LUT_B34(xx(:));
            Bm34(:)             = LUT_Bm34(xx(:));

            %f_0 (Dinardo et al. (2018), Eq. 23)
            % f0   = pi/(2*sqrt(2))*abs((t-X2)./(2*sigma_c)).^0.5.*...
            %     (besseli(-1/4,((t-X2)./(2*sigma_c)).^2,1)+...
            %     (t-X2)./abs(t-X2).*besseli(1/4,((t-X2)./(2*sigma_c)).^2,1));
            f0_0 = (pi*2^(3/4)) / (4 * gamma(3/4));
            f0   = pi/(2*sqrt(2))*abs((t-X2)./(2*sigma_c)).^0.5.*...
                (Bm14 + (t-X2)./abs(t-X2).*B14);
            f0(xx < eps) = f0_0;

            %f_1 (Dinardo et al. (2018), Eq. 24)
            % f1   = pi/(2*sqrt(2))*abs((t-X2)./(2*sigma_c)).^1.5.*...
            %     ((besseli(1/4,((t-X2)./(2*sigma_c)).^2,1)-...
            %     besseli(-3/4,((t-X2)./(2*sigma_c)).^2,1))+...
            %     (t-X2)./abs(t-X2).*(besseli(-1/4,((t-X2)./(2*sigma_c)).^2,1)...
            %     -besseli(3/4,((t-X2)./(2*sigma_c)).^2,1)));
            f1_0 = -(2^(3/4) * gamma(3/4)) / 4;
            f1   = pi/(2*sqrt(2))*abs((t-X2)./(2*sigma_c)).^1.5.*...
                ((B14-Bm34) + (t-X2)./abs(t-X2).*(Bm14-B34));
            f1(xx < eps) = f1_0;
            
            %Single look waveform (Dinardo et al. (2018), Eq. 22)
            W_sl = (X1*sqrt(2*pi)*DDAcf.alpha_PF^2)*gamma_0.*(f0+SmodID*sigma_z/L_gamma*T_k.*g_l*sigma_s.*f1)./(DDAcf.B*sigma_c).^0.5;

            %Modelled multi-looked SAR (Delay-Doppler) echo (Dinardo et al. (2018), Eq. 28)
            y = TN + mean(W_sl,2);
        end
        
        function y = SAMOSA_DPM_v2_5_2_fun(x,n,DDAcf,LUT_AP,LUT_F0,LUT_F1,BeamIDX)
            %3/4-parameter SAMOSA2/SAMOSA3 retrackers described in
            %Gommenginger et al. 2017. Sometimes, the version of the DPM is
            %explicitly mentioned.
            %
            %Gommenginger, C., Martin-Puig, C., Srokosz, M., Caparrini, M.,
            %Dinardo, S., Lucas, B., Restano, M., Ambrozio, A. & J.
            %Benveniste (2017). Detailed Processing Model of the Sentinel-3
            %SRAL SAR altimeter ocean waveform retracker, Version 2.5.2, 31
            %October 2017. ESA-ESRIN Contract No. 20698/07/I-LG
            %"Development of SAR Altimetry Mode Studies and Applications
            %over Ocean, Coastal Zones and Inland Water" (SAMOSA), 85
            %pages.
            %
            %Input: 
            % -x        = vector with unknowns;
            % -n        = bin indices WF, TN, SWH, loc. radius of curvature 
            %             of the Earth's surface, dist. satellite-surface,
            %             tangential velocity, slope of orbit, pitch & roll
            %             angles, SmodID, and EstPar
            % -DDAcf    = DDA configuration file
            % -LUT_AP   = look-up-table alph_P
            % -LUT_F0   = look-up-table f0 term
            % -LUT_F1   = look-up-table f1 term
            % -BeamIDX  = vector with Doppler beam indices
            %
            % Output:
            % -y      = the SAMOSA(2/3)(FF) model
            %
            % Author: D.C. Slobbe (February 2019)
            
            %Unpack vector n
            EstPar    = n(end);
            SmodID    = n(end-1);
            ksiy      = n(end-2); %roll mispointing angle (across-track)
            ksix      = n(end-3); %pitch mispointing angle (along-track)
            Orb_slope = n(end-4);
            vs        = n(end-5);
            h         = n(end-6);
            Re        = n(end-7);
            DUMx      = n(end-8);
            TN        = n(end-9);
            n         = n(1:end-10);

            %Unpack vector with unknowns. Convert X2 to seconds
            if EstPar == 3
                X1 = x(1); X2=x(2)/1E9; X3=x(3); X4=0;
            elseif EstPar == 4
                X1 = x(1); X2=x(2)/1E9; X3=DUMx; X4=x(end)*1E5;
            else
                X1 = x(1); X2=x(2)/1E9; X3=x(3); X4=x(4)*1E5;
            end

            %Preliminaries
            Ns        = length(n);                   %Nr of bins/samples in any waveform (after zero-padding)
            os_ZP     = Ns/DDAcf.Np;                 %Zero-Padding Oversampling Factor
            RefBin    = (DDAcf.RefBin-1)*os_ZP + 1;  %Reference bin

            %SL2. Compute the spherical surface parameter
            alpha     = 1+h/Re;                      %Eq. 58

            %ML4. Sub-setting the Stack
            dfa                = DDAcf.fp/DDAcf.Nb;                       %Eq. 61
            Dopp_Freq_Stack_TS = BeamIDX*dfa;
            Neff               = numel(Dopp_Freq_Stack_TS);
            fa                 = Dopp_Freq_Stack_TS;

            %SL3. Calculate some new parameters
            dtau      = 1/DDAcf.B;                                %Time resolution [s] (Eq. 60)
            t_s       = dtau/os_ZP;                               %Waveform sampling interval [s] (~= dtau because of zero-padding)
            xp        = h*tan(ksix);                              %Eq. 62
            yp        = -h*tan(ksiy);                             %Eq. 63
            alphax    = DDAcf.shx*8*log(2)/(h^2*DDAcf.theta_x^2); %Eq. 64
            alphay    = DDAcf.shy*8*log(2)/(h^2*DDAcf.theta_y^2); %Eq. 65
            Tb        = DDAcf.Nb/DDAcf.fp;                        %Eq. 66
            Lx        = (DDAcf.c*h)/(2*vs*DDAcf.fc*Tb);           %Eq. 67
            Ly        = sqrt((DDAcf.c*h)/(alpha*DDAcf.B));        %Eq. 68
            Lz        = DDAcf.c/(2*DDAcf.B);                      %Eq. 69
            Lg        = alpha/(2*h*alphay);                       %Eq. 70
            sigmaz    = X3/4;                                     %Eq. 71
 
            %SL4. Define indices in delay (K) and doppler (L) space
            t         = n*t_s-t_s*RefBin;
            tau       = t-X2;                        %Eq. 72
            K         = tau/dtau;                    %Eq. 73
            L         = fa/dfa;                      %Eq. 74
 
            %SL5. Define a value for alph_P
%             if EstPar == 4
%                 alph_P = 0.552;
%             else
            alph_P    = LUT_AP(X3);
%             end

            %SL5. Calculate GL and GLK
            Ls    = Orb_slope * h / (alpha * Lx);
            gamma = 2 * (L-Ls) * Lx^2 / Ly^2;        %Eq. 75
            if X3 ~= 0
                sign = X3/abs(X3);
            else
                sign = 1;
            end
%                 sigma_p = alph_P/DDAcf.B;
%                 theta_l = Lx/h*(L-Ls);
%                 GL  = 1./(DDAcf.B*sqrt(sigma_p^2*(1+(Lx/Lz)^2*(alpha*theta_l).^2)+(2*sigmaz/DDAcf.c)^2));
            GL  = 1 ./ sqrt( (alph_P)^2 + (alph_P)^2*gamma.^2 + sign*(X3/(4*Lz))^2 ); %Eq. 78
            GLK = GL .* K;                                                            %Eq. 79

            %SL6. Calculate GAMMA0 and Const
            XL         = L*Lx;                       %Eq. 80
            YK         = Ly*sqrt(K);                 %Eq. 82
            YK(K <= 0) = 0;                          %Eq. 81
            GAMMA0     = exp( -alphay*yp^2 - alphax*(XL-xp).^2 - XL.^2*X4/h^2 ) .* exp( -(alphay+X4/h^2)*YK.^2) .* cosh(2*alphay*yp*YK); %Eq. 83
            Const      = (DDAcf.alpha_PF)^2*sqrt(2*pi);

            %SL6. Calculate Tk
            Tk         = (1 + X4/(alphay*h^2)) - (yp/(Ly*sqrt(K)))'.*tanh(2*alphay*yp*Ly*sqrt(K));
            Tk(K <= 0) = (1 + X4/(alphay*h^2)) - 2*alphay*yp^2;

            %Obtain LUT_F0(GLK) and LUT_F1(GLK)
            [F0,F1]          = deal(nan(Ns,Neff));
            F0(:)            = LUT_F0(GLK(:));
            F0(GLK < -19)    = 0;
            F0(GLK > 41.999) = 1/2*sqrt(2*pi./GLK(GLK > 41.999));
            F1(:)            = LUT_F1(GLK(:));
            F1(GLK < -19)    = 0;
            F1(GLK > 41.999) = 1/4*sqrt(2*pi./GLK(GLK > 41.999).^3);
            
            %SL7. Calculate the final expression of Pr_SL
            Pr_Stack     = X1 * Const .* sqrt(GL) .* GAMMA0 .*(F0+SmodID*(sigmaz/Lg)*(sigmaz/Lz)*GL.*Tk.*F1); %Eq. 84
            
            %ML6. Mask out from the Stack of Doppler beams the power bins located beyond the radar window length
            Lx           = (DDAcf.lambda0 * h) / (2 * vs * Tb);
            Beam_Range   = repmat(h*( sqrt( 1+alpha*(Lx*BeamIDX/h).^2 ) -1 ),Ns,1);
            dr           = repmat(((Ns-1):-1:0)',1,Neff);
            Window_Range = DDAcf.c/(2*DDAcf.B)*dr;
            Pr_Stack(Beam_Range > Window_Range) = 0;
            
            %ML7. Apply Weighting Function to the Stack prior to incoherent summation
            %PR_Stack = Pr_Stack .* Stack_Weights(j); %Eq. 91
 
            %ML8. Calculate multi-looked waveform by incoherent summation across the Stack
            y = TN + sum(Pr_Stack,2); %Eq. 92
        end
        
        function [n_ret,A] = Threshold(WD,sub_n,sub_WD)
            %Threshold retracker designed by Davis (Davis, C. H. (1997). A
            %robust threshold retracking algorithm for measuring ice-sheet
            %surface elevation change from satellite radar altimeters.
            %Geoscience and Remote Sensing, IEEE Transactions on, 35(4),
            %974-979.)
            %
            %Input: 
            % -WD     = full waveform;
            % -sub_n  = bin indices of sub-waveform;
            % -sub_WD = sub-waveform.
            %
            % Output:
            % -n_ret  = retracking location [bins].
            %
            % Author: D.C. Slobbe, January 2016

            %Settings
            nalias = 12;  %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)
            Tcoeff = 0.6; %Percentage of the max waveform amplitude above DC/noise level
            
            %Maximum amplitude of sub-waveform
            Amax   = max(sub_WD);
            
            %Compute thermal noise level
            DC     = mean(WD(1:nalias));

            %Compute threshold level
            TL     = DC + Tcoeff*(Amax-DC);
            
            %Compute location of first gate that exceeds threshold level
            [~,ID] = max(sub_WD>TL,[],1);
            n_hat  = sub_n(sub2ind(size(sub_n),1:size(sub_n,1),ID));
            
            %Compute precise retracking location on the leading edge of the
            %waveform by linear interpolation between the bins adjacent to
            %the threshold crossing
            n_ret = (n_hat - 1) + (TL - WD(n_hat - 1)')./(WD(n_hat) - WD(n_hat - 1))';
            n_ret(WD(n_hat) - WD(n_hat - 1) <= 2*DC) = n_hat(WD(n_hat) - WD(n_hat - 1) <= 2*DC) - 1;
            
            %Amplitude (needed to compute sigma0)
            A     = Amax-DC;
        end
        
        function [n_ret,A] = Threshold_mat(WD)
            %Threshold retracker designed by Davis (Davis, C. H. (1997). A
            %robust threshold retracking algorithm for measuring ice-sheet
            %surface elevation change from satellite radar altimeters.
            %Geoscience and Remote Sensing, IEEE Transactions on, 35(4),
            %974-979.)
            %
            %Input: 
            % -WD     = set of full waveforms or sub-waveforms;
            %
            % Output:
            % -n_ret  = retracking location [bins].
            %
            % Author: D.C. Slobbe, January 2019

            %Settings
            nalias = 12;  %Number of unaliased waveform samples used to estimate DC level (i.e. thermal noise)
            Tcoeff = 0.6; %Percentage of the max waveform amplitude above DC/noise level
            
            %Maximum amplitude of sub-waveform
            Amax   = max(WD);
            
            %Compute thermal noise level
            DC     = mean(WD(1:nalias,:));

            %Compute threshold level
            TL     = DC + Tcoeff*(Amax-DC);
            
            %Compute location of first gate that exceeds threshold level
            [~,ID] = max(WD>TL,[],1);
            n_hat  = ID;
            n_ID   = sub2ind(size(WD),ID,1:size(WD,2));
            
            %Compute precise retracking location on the leading edge of the
            %waveform by linear interpolation between the bins adjacent to
            %the threshold crossing
            n_ret = (n_hat - 1) + (TL - WD(n_ID - 1))./(WD(n_ID) - WD(n_ID - 1));
            n_ret(WD(n_ID) - WD(n_ID - 1) <= 2*DC) = n_hat(WD(n_ID) - WD(n_ID - 1) <= 2*DC) - 1;
            
            %Amplitude (needed to compute sigma0)
            A     = Amax-DC;
        end
        
    end
end
