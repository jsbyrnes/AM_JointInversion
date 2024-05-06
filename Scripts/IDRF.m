function [RF, RF_Time] = IDRF(Phase,P,D,dt,t_for,t_trun,gauss_t,accept_mis1,accept_mis2,itmax)
% Iterative Deconvolution and Receiver-Function Estimation in time domain
% Version used in Junlin Hua's Sp RF + CCP stacking code, with some cleaning.

% Phase is 'Sp' or 'Ps'
% P for parent phase (S in Sp case), D for daughter phase.
% dt: time between samples
% gauss_t: 1 std for the Gaussian convolved
% accept_mis1, accept_mis2 and itmax for stop criteria of the loop
% t_for,t_trun: remove amplitude with t > 0 in the Sp case 

    misfit = 1;
    misfit_old = 9999999999999;
    misfit_ref = sqrt(sum(D.^2));

    RF_tmp = zeros(length(P)*2-1,1);

    D_cur = D;

    itnum = 0;

    [~,t_corr] = xcorr(D_cur,P);
    RF_Time = t_corr*dt;

    while (misfit_old-misfit>accept_mis1*misfit || misfit>accept_mis2*misfit_ref) && itnum <= itmax
        [amp_corr,~] = xcorr(D_cur,P);
        auto_corr = xcorr(P);
        [~,ind] = max(abs(amp_corr));
        amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
        
        RF_tmp(ind) = RF_tmp(ind)+amp_rf;
        D_sub = conv(P,RF_tmp,'same');
        D_cur = D - D_sub;
        
        misfit_old = misfit;
        misfit = sqrt(sum(D_cur.^2))/misfit_ref;
        itnum = itnum+1;
    end

    RF = RF_tmp;

    if strcmp(Phase,'Sp')
        RF(RF_Time>t_trun)=0;
        RF = RF(RF_Time<=t_for);
        RF_Time = RF_Time(RF_Time<=t_for);
    else
        RF(RF_Time<t_trun)=0;
        RF = RF(RF_Time>=t_for);
        RF_Time = RF_Time(RF_Time>=t_for);
    end

    if gauss_t~=0
%         gauss_sig = gauss_t/dt;
%         x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
%         Gauss_win = exp(-x.^2/(2*gauss_sig^2));
%         %Gauss_win = exp(-abs(x)/gauss_sig);
% 
%         %C = exp_cov(length(RF), gauss_t/dt);
%         %RF = chol(C)*(RF.*tukeywin(length(RF), 0.1));
% 
%         RF = conv(RF,Gauss_win,'same');

%         Nfft                   = length(RF);         % number  of points in fft = n of samples
%         dF                     = 1/(dt*Nfft);       % frequency interval
%         Nyq                    = (1/(2*dt));        % Nyquist frequency
%         freqVals               = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
%         freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);
%         %freqVals               = abs(freqVals);
%         
%         %Faux = gauss_t;%number is half width in s
%         A    = exp(-((freqVals.^2)/(2*gauss_t^2)));
%         %A    = exp(-((abs(freqVals - gauss_t))/(Faux)));
%         A    = A'; %for dimensional consistency with inTr.data
%         
%         %Apply to the traces and revert to time domain
%         
%         RF  = real(ifft(fft(RF).*A));
% 
%         %restore the amplitude
% 
%         impulse = zeros(size(RF));
%         impulse(round(length(impulse)/2)) = 1;
%         impulse  = real(ifft(fft(impulse).*A));
% 
%         RF = RF/max(impulse);

    end
    
end

