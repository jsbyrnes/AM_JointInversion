function [ trace, trace1 ] = multitaper2rf( P_comp, SV_comp, dt, rf_shift, ntap, NW, nmtw, phase, filter_limits )
%MULTITAPER2RF Use with synthetic RF traces to make quick receiver
%functions striaght from the synthetic source
%byrnes.joseph@gmail.com

    [E,V]=dpss(ntap, NW, nmtw);
    
    C.rf_shift = rf_shift;
    C.ntap = ntap;
    C.nmtw = nmtw;
        
    if phase == 'P'
        
        G = drive_mtm(P_comp, SV_comp, zeros(size(SV_comp)), C, E, V, 1/dt);
        
    elseif strcmp(phase, 'SV')
        
        %need to reverse the time, clip, and flip P polarity
        
        %first get the main arrival time and length, both in samples
        nsamples = length(SV_comp);
        
        SV_comp = flipud(SV_comp);
        P_comp = flipud(-1*P_comp);

        [~, arrival_index] = max(SV_comp);
        
        SV_comp = SV_comp(arrival_index - C.rf_shift/dt + 1:end);
        P_comp = P_comp(arrival_index - C.rf_shift/dt + 1:end);
        
        %how many zeros are needed to get back original length?
        new_nsamples = length(SV_comp);
        
        SV_comp = [SV_comp; zeros(nsamples - new_nsamples, 1)];
        P_comp = [P_comp; zeros(nsamples - new_nsamples, 1)];
        
        G = drive_mtm(SV_comp, P_comp, zeros(size(P_comp)), C, E, V, 1/dt);
    
    end
    
    %detrend the RFs just in case
    
    G.C1 = detrend(G.C1);
    G.C2 = detrend(G.C2);
    G.C3 = detrend(G.C3);
            
    trace1 = bandpassfilt_rfs(G.C1, dt, filter_limits(1), filter_limits(2));
    trace = bandpassfilt_rfs(G.C2, dt, filter_limits(1), filter_limits(2));
    
    %normalize to seis1
            
    index = rf_shift/dt;
    
    norm_window = index + [-3 3]/dt;
    
    amp = max(trace1(norm_window(1):norm_window(2)));
    
    trace1 = trace1/amp;
    trace = trace/amp;
        
end

