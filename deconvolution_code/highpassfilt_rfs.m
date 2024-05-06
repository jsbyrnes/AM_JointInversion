function [ ts_out ] = highpassfilt_rfs( ts, dt, corner, order )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a calculated r value

    n = length(ts);
    
    %tap = tukeywin(n, 1/(n*dt*low));
    tap = tukeywin(n, 0.1);
    
    if nargin == 5
        
        [z,p,k] = butter(order,corner.*(2*dt), 'high');
        
    else
        
        [z,p,k] = butter(4 ,corner.*(2*dt), 'high');

    end
    
    [SOS, G] = zp2sos(z,p,k);
    
    ts_out = filtfilt(SOS, G, ts.*tap);
    
end

