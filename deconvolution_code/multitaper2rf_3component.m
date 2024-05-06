function [ trace_p, trace_sv, trace_sh, sv_qual] = multitaper2rf_3component( P_comp, SV_comp, SH_comp, ...
    dt, rf_shift, E, V, phase, filter_limits, taper_ind, arrival_time)
%MULTITAPER2RF Use with synthetic RF traces to make quick receiver
%functions striaght from the synthetic source. Set up for three component
%recevier functions with P
%byrnes.joseph@gmail.com

    %[E,V]=dpss(ntap, NW, nmtw);
    
    old_n = length(P_comp);

%     tap = tukeywin(old_n, 0.25);
% 
%     P_comp  = [ zeros(size(P_comp));  P_comp.*tap;   zeros(size(P_comp))];
%     SV_comp = [ zeros(size(SV_comp)); SV_comp.*tap;  zeros(size(SV_comp))];
%     SH_comp = [ zeros(size(SH_comp)); SH_comp.*tap;  zeros(size(SH_comp))];

    P_comp  = [ zeros(size(P_comp));  P_comp;   zeros(size(P_comp))];
    SV_comp = [ zeros(size(SV_comp)); SV_comp;  zeros(size(SV_comp))];
    SH_comp = [ zeros(size(SH_comp)); SH_comp;  zeros(size(SH_comp))];

    C.rf_shift = round(dt*length(P_comp)/2);
    [ C.ntap, C.nmtw ] = size(E);
        
    if phase == 'P'
        
        G = drive_mtm(P_comp, SV_comp, SH_comp, C, E, V, 1/dt);
        
        sv_reduction = (SV_comp - conv(P_comp, G.C2, 'same')).^2;
        sv_reduction = sum(sv_reduction(old_n:2*old_n));

        SV_power = SV_comp.^2;
        SV_power = sum(SV_power(old_n:2*old_n));

        sv_qual = 1 - sv_reduction/( SV_power );

        %detrend the RFs just in case

%         G.C1 = detrend(G.C1);
%         G.C2 = detrend(G.C2);
%         G.C3 = detrend(G.C3);

        %trace_p = bandpassfilt_rfs(G.C1, dt, filter_limits(2), filter_limits(1));
        %trace_sv = bandpassfilt_rfs(G.C2, dt, filter_limits(2), filter_limits(1));
        %trace_sh = bandpassfilt_rfs(G.C3, dt, filter_limits(2), filter_limits(1));

        trace_p  = lowpassfilt_rfs(G.C1, dt, filter_limits(2));
        trace_sv = lowpassfilt_rfs(G.C2, dt, filter_limits(2));
        trace_sh = lowpassfilt_rfs(G.C3, dt, filter_limits(2));

        %normalize to seis1

        index = round(C.rf_shift/dt);

        norm_window = round(index + [-rf_shift/2 rf_shift/2]/dt);

        amp = max(trace_p(norm_window(1):norm_window(2)));

        trace_p = trace_p/amp;
        trace_sv = trace_sv/amp;
        trace_sh = trace_sh/amp;
        
        %and now shift back for display
        trace_p  = circshift(trace_p, round(rf_shift/dt) - index);
        trace_sv = circshift(trace_sv, round(rf_shift/dt) - index);
        trace_sh = circshift(trace_sh, round(rf_shift/dt) - index);

        trace_p  = trace_p (1:old_n);%.*tukeywin(old_n, 0.1);
        trace_sv = trace_sv(1:old_n);%.*tukeywin(old_n, 0.1);
        trace_sh = trace_sh(1:old_n);%.*tukeywin(old_n, 0.1);

    elseif strcmp(phase, 'SV') || strcmp(phase, 'S')
        
        %need to reverse the time, clip, and flip P polarity
        
        %first get the main arrival time and length, both in samples
        nsamples = length(SV_comp);
        
        %taper the SV comp
        SV_comp = taper_ts(SV_comp, taper_ind(1), taper_ind(2), 0.3);        
        
        SV_comp = flipud(SV_comp);
        P_comp = flipud(-1*P_comp);

        if ~isempty(arrival_time)
        
            arrival_index = round(arrival_time/dt);%preflipping
            arrival_index = length(SV_comp) - arrival_index;
            
        else
            
            [~, arrival_index] = max(abs(SV_comp));%only appropriate if you used tapering
            
        end
        
        try
        
            SV_comp = SV_comp(arrival_index - round(C.rf_shift/dt) + 1:end);
            P_comp = P_comp(arrival_index - round(C.rf_shift/dt) + 1:end);

        catch
            
            SV_comp = SV_comp(round(C.rf_shift/dt) + 1:end);
            P_comp = P_comp(round(C.rf_shift/dt) + 1:end);
            
        end
        
        %how many zeros are needed to get back original length?
        new_nsamples = length(SV_comp);
        
        SV_comp = [SV_comp; zeros(nsamples - new_nsamples, 1)];
        P_comp = [P_comp; zeros(nsamples - new_nsamples, 1)];
        
        SV_comp = SV_comp - mean(SV_comp);
        P_comp = P_comp - mean(P_comp);
        
        G = drive_mtm(SV_comp, P_comp, zeros(size(P_comp)), C, E, V, 1/dt);
    
        %detrend the RFs just in case

        G.C1 = detrend(G.C1);
        G.C2 = detrend(G.C2);
        G.C3 = detrend(G.C3);

        trace_sv = bandpassfilt_rfs(G.C1, dt, filter_limits(2), filter_limits(1));
        trace_p = bandpassfilt_rfs(G.C2, dt, filter_limits(2), filter_limits(1));
        trace_sh = bandpassfilt_rfs(G.C3, dt, filter_limits(2), filter_limits(1));

        %normalize to seis1

        index = round(rf_shift/dt);

        norm_window = round(index + [-rf_shift/2 rf_shift/2]/dt);

        amp = max(abs(trace_sv(norm_window(1):norm_window(2))));

        trace_sv = trace_sv/amp;
        trace_p = trace_p/amp;
        trace_sh = trace_sh/amp;

    elseif strcmp(phase, 'SH')
        
        %need to reverse the time and clip
        
        %first get the main arrival time and length, both in samples
        nsamples = length(SH_comp);
        
        %taper the SH comp
        SH_comp = taper_ts(SH_comp, taper_ind(1), taper_ind(2), 0.3);
        
        SH_comp = flipud(SH_comp);
        P_comp = flipud(P_comp);
        
        [~, arrival_index] = max(SH_comp);
        
        try
        
            SH_comp = SH_comp(arrival_index - round(C.rf_shift/dt) + 1:end);
        
        catch
            
            keyboard
            
        end
        
        P_comp = P_comp(arrival_index - round(C.rf_shift/dt) + 1:end);
        
        %how many zeros are needed to get back original length?
        new_nsamples = length(SH_comp);
        
        SH_comp = [SH_comp; zeros(nsamples - new_nsamples, 1)];
        P_comp = [P_comp; zeros(nsamples - new_nsamples, 1)];
        
        G = drive_mtm(SH_comp, P_comp, zeros(size(P_comp)), C, E, V, 1/dt);

        %detrend the RFs just in case

        G.C1 = detrend(G.C1);
        G.C2 = detrend(G.C2);
        G.C3 = detrend(G.C3);

        trace_sh = bandpassfilt_rfs(G.C1, dt, filter_limits(2), filter_limits(1));
        trace_p = bandpassfilt_rfs(G.C2, dt, filter_limits(2), filter_limits(1));
        trace_sv = bandpassfilt_rfs(G.C3, dt, filter_limits(2), filter_limits(1));

        %normalize to seis1

        index = round(rf_shift/dt);

        norm_window = round(index + [-rf_shift/2 rf_shift/2]/dt);

        amp = max(trace_sh(norm_window(1):norm_window(2)));

        trace_sh = trace_sh/amp;
        trace_p = trace_p/amp;
        trace_sv = trace_sv/amp;
        
    end
            
end

