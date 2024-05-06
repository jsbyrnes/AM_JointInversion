function [ ts ] = bandpassfilt( ts, dt, high, low )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a r = .2 Tukey window(hardwired, but easy to
%change)

    n = length(ts);
    
    tap = tukeywin(n, .2);
    
    [z,p,k] = butter(3,[low high].*(2*dt));

    [SOS,G] = zp2sos(z,p,k);
    
    ts = filtfilt(SOS,G,ts.*tap);
    
end

