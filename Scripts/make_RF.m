function rf_out = make_RF(P, D, dt, gsig, pre, filterstyle)

    %pad
    P = [ zeros(length(P), 1); P; zeros(length(P), 1) ];
    D = [ zeros(length(D), 1); D; zeros(length(D), 1) ];

    if strcmp(filterstyle, 'gaussian')

        gsig = 1/(pi*gsig);

        gauss_sig = round(gsig/dt);
    
        x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
        Gauss_win = exp(-x.^2/(2*gauss_sig^2));
        P = conv(P,Gauss_win,'same');
        D = conv(D,Gauss_win,'same');

    elseif strcmp(filterstyle, 'butter')

        P = bandpassfilt_rfs(P, dt, gsig, 0.01);
        D = bandpassfilt_rfs(D, dt, gsig, 0.01);                

    end

    %hardwired, oh well
    [rf_out, ~] = IDRF('P', P, D, dt, -pre, -pre, gsig, 0.001, 0.001, 400);
        
    if strcmp(filterstyle, 'gaussian')
    
        rf_out = conv(rf_out,Gauss_win,'same');
        %D = conv(D,Gauss_win,'same');

    elseif strcmp(filterstyle, 'butter')

        rf_out = bandpassfilt_rfs(rf_out, dt, gsig, 0.01);

        %impulse for normalization
        impulse = zeros(size(rf_out));
        impulse(round(length(impulse)/2)) = 1;
        impulse = bandpassfilt_rfs(impulse, dt, gsig, 0.01);

        rf_out = rf_out/max(impulse);

    end

end


    %gsig = gsig/(2*pi);
    %gsig = gsig/sqrt(2);

%     Nfft                   = length(P);         % number  of points in fft = n of samples
%     dF                     = 1/(dt*Nfft);       % frequency interval
%     Nyq                    = (1/(2*dt));        % Nyquist frequency
%     freqVals               = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
%     freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);
%     %freqVals               = abs(freqVals);
%     
%     A    = exp(-(freqVals.^2)/(2*gsig^2));
%     %A    = A'; %for dimensional consistency with inTr.data
%     
%     %Apply to the traces and revert to time domain
%     P  = real(ifft(fft(P).*A));
%     D  = real(ifft(fft(D).*A));