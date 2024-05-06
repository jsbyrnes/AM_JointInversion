function [ ts_out ] = bandpassfilt_rfs( ts, dt, high, low )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a calculated r value

    n = length(ts);
    
    tap = tukeywin(n, 1/(n*dt*low));
    
    [z,p,k] = butter(3,[low high].*(2*dt));
    %[z,p,k] = ellip(3, 0.001, 100, [low high].*(2*dt));
    %[z,p,k] = cheby2(3, 50, [low high].*(2*dt));
    
    [SOS, G] = zp2sos(z,p,k);
    
    ts_out = filtfilt(SOS, G, ts.*tap);

end

