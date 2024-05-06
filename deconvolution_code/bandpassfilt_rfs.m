function [ ts_out ] = bandpassfilt_rfs( ts, dt, high, low, order )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a calculated r value

    n = length(ts);
    
    %tap = tukeywin(n, 1/(n*dt*low));
    tap = tukeywin(n, 0.33);
    
    if nargin == 5
        
        [z,p,k] = butter(order,[low high].*(2*dt));
        
    else
        
        [z,p,k] = butter(4 ,[low high].*(2*dt));

    end
    
    [SOS, G] = zp2sos(z,p,k);
    
    ts_out = filtfilt(SOS, G, ts.*tap).*tap;
    
end

