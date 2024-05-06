function [ ts_out ] = cheby1filt_rfs( ts, dt, high, low, order, ripple_db )
%BANDPASSFILT Bandpass filters the data between high and low in Hz
%Data is tapered first with a calculated r value

    n = length(ts);
    
    %tap = tukeywin(n, 1/(n*dt*low));
    tap = tukeywin(n, 0.1);
            
    %[z,p,k] = ellip(order, ripple_db, stop_attentuation, [low high].*(2*dt));
    [z,p,k] = cheby1(order, ripple_db, [low high].*(2*dt));
    
    [SOS, G] = zp2sos(z,p,k);
    
    ts_out = filtfilt(SOS, G, ts.*tap);
    
end

