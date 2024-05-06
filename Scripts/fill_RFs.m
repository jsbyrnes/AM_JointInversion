function allWfs = fill_RFs(vs_surface, vp_surface, t_interp, rawData, Parameters, kbts)
    %The window structure here is very hardwired. Oh well. 

    %%%%%%%%%%%%
    %create the allWfs struct
    for i1 = 1:length(Parameters.p_array)
    
        for i2 = 1:length(Parameters.g_array)
    
            allWfs(i1, i2).p   = Parameters.p_array(i1);
            allWfs(i1, i2).g   = Parameters.g_array(i2);
            allWfs(i1, i2).n   = 0;
            allWfs(i1, i2).rfr = [];
            allWfs(i1, i2).rft = [];
            allWfs(i1, i2).Z   = [];
            allWfs(i1, i2).R   = [];
    
        end
    
    end

    G     = length(Parameters.g_array);

    weights = zeros(length(Parameters.p_array), G);

    for j = 1:G

        window = max([ 3200 - 800*(j-1), 800]);
        
        %window = 400;
        [E,V] = dpss(window/8, 2.5, 3);

        for k = 1:length(rawData)            
    
            if ~isnan(vs_surface)

                [ P_comp, SV_comp ] = ZR2PSV(rawData(k).Z(1:window), -1*rawData(k).R(1:window), ...
                    vp_surface, vs_surface, rawData(k).p);%surface vp won't matter
                P_comp  = -1*P_comp;
                SV_comp = -1*SV_comp;

            else

                P_comp  = rawData(k).Z(1:window);
                SV_comp = rawData(k).R(1:window);

            end
       
            SH_comp = zeros(size(SV_comp));%need to flip?;
                            
            filterlims = [ 0.05 Parameters.g_array(j) ];
    
            if strcmp(Parameters.rf_style, 'multi-taper')

                try

                    [~, rfr, ~, ~] = multitaper2rf_3component(P_comp, SV_comp, ...
                        SH_comp, 1/Parameters.sample_rate, abs(Parameters.datawin(1)), E, V, 'P', ...
                        filterlims, [ 1 length(P_comp)]);

                catch

                    keyboard

                end

            elseif strcmp(Parameters.rf_style, 'time-domain')

                rfr = make_RF(P_comp, SV_comp, 0.05, Parameters.g_array(j), abs(Parameters.datawin(1)), Parameters.filter_style);

            end
                
            t   = ((0:length(rfr)-1))*0.05 + Parameters.datawin(1);
            
            rfr = interp1(t, rfr, t_interp, 'linear');%force align but not a resample
            
            if isrow(rfr)
            
                rfr = rfr';
            
            end
                                
            %here the clip to what we will invert

            ind = t_interp >= -Parameters.pre & t_interp <= Parameters.max_time;

            rawData(k).rfr(:, j)   = rfr(ind).*tukeywin(sum(ind), 0.1);
        
        end

    end

    Parameters.t = ((-Parameters.pre:0.05:Parameters.max_time))';

    %stack and normalize
    for j = 1:G
        
        for tix = 1:length(Parameters.t)

            for k = 1:length(rawData)
    
                rfr = rawData(k).rfr(:, j);
                dp  = rawData(k).p - Parameters.move_out.ref_p;

                %and reinterpolate to new time
                rfr          = interp1(Parameters.t*(1 - dp*Parameters.move_out.time_slope), rfr, Parameters.t, 'spline', 'extrap');

                rfr = rfr*(1 - dp*Parameters.move_out.amp_slope);
                
                rfr_val(k) = rfr(tix);
                baz_val(k) = rawData(k).baz;
                snr_val(k) = rawData(k).rqual;

            end

            if isempty(Parameters.baz)

                opt = optimoptions(@fminunc, 'Display', 'none');
    
                parfor k = 1:kbts
    
                    if k == 1
    
                        ind = 1:length(rfr_val);
    
                    else
    
                        ind = randi(length(rfr_val), [ length(rfr_val) 1]);
    
                    end
    
                    if isempty(Parameters.baz_range)
    
                        cosfit = @(a) sum( (((rfr_val(ind) - (a(1) + a(2)*cosd(2*baz_val(ind)) + a(3)*sind(2*baz_val(ind)) + a(4)*cosd(baz_val(ind)) + a(5)*sind(baz_val(ind))))/0.05).^2).*sqrt(snr_val))/sum(sqrt(snr_val)) ...
                            + (a(2)^2 + a(3)^2)/(0.05 + (length(rfr_val)^2)/2e4)^2 + (a(4)^2 + a(5)^2)/(0.05 + (length(rfr_val)^2)/2e4)^2;%damp if there aren't many. Weighting does very little. 
        
                        x = fminunc(cosfit, [mean(rfr_val(ind));0;0;0;0], opt);
            
                        val(k)     = x(1);
                        val_amp1(k) = sqrt(x(4)^2 + x(5)^2);
                        val_amp2(k) = sqrt(x(2)^2 + x(3)^2);
    
                    else
    
                        val(k) = mean(rfr_val(ind));
    
                    end
    
                end

                allWfs(1, j).rfr(tix, 1)         = mean(val);
                allWfs(1, j).rfr_1T(tix, 1)      = mean(val_amp1);
                allWfs(1, j).rfr_2T(tix, 1)      = mean(val_amp2);
                allWfs(1, j).rfr_std(tix, 1)     = std(val);% + Parameters.hd_error*(allWfs(1, j).rfr_1T(tix, 1) + allWfs(1, j).rfr_2T(tix, 1));
                %allWfs(1, j).rfr_1T_std(tix, 1)  = std(val_amp1);
                allWfs(1, j).rfr_2T_std(tix, 1)  = std(val_amp2);

            else

                allWfs(1, j).rfr(tix, 1)         = mean(rfr_val);
                allWfs(1, j).rfr_std(tix, 1)     = std(rfr_val)/sqrt(length(rfr_val));

            end

            allWfs(1,j).n                    = 1;

        end

        if isempty(Parameters.baz)

            aniso = smooth(allWfs(1,j).rfr_2T + allWfs(1,j).rfr_1T, 20/Parameters.g_array(j));
            aniso = aniso - 0.03;
            aniso(aniso<0) = 0;
            allWfs(1,j).rfr_std = allWfs(1,j).rfr_std + 0.03*(1 - exp(-aniso.^2/(0.5*0.03^2)));
            
        end
        
        allWfs(1,j).rfr_std( allWfs(1,j).rfr_std < 0.01*Parameters.min_error*max(abs(allWfs(1, j).rfr))) = 0.01*Parameters.min_error*max(abs(allWfs(1, j).rfr));
        
        %and again, set a floor. 
        allWfs(1,j).rfr_std( allWfs(1,j).rfr_std < 0.005) = 0.005;

    end
    
end

            %gauss_sig = (1/Parameters.g_array(j))/0.05;
            %x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
            %Gauss_win = exp(-x.^2/(2*gauss_sig^2));
            %P_comp    = conv(P_comp,Gauss_win,'same');
            %SV_comp   = conv(SV_comp,Gauss_win,'same');
            
            %[~, ind] = max(P);
            
            %rf = rf/max(abs(P));
            %P  = P/max(abs(P));
            
            %filter_limits = [ 0.01 Parameters.g_array(j)/pi ];
            
            %P  = bandpassfilt_rfs(P_comp + randn(size(P_comp))*noiseamp,  0.05, filter_limits(2), filter_limits(1));
            %SV = bandpassfilt_rfs(SV_comp + randn(size(P_comp))*noiseamp, 0.05, filter_limits(2), filter_limits(1));

%     %%%%%%%%%%%%%%%%%%%%%
%     %check if fader helps
%     Conclusion - on synthetic, it does an excellent job of removing the
%     reverbs and isolating the moho. But, then there is a giant impulse at
%     short lag! how am I suppose to invert that?

%     [SedPre,tlaga,r0,tac,ac,dsin] = DeterminationTest_func(fader_array,Parameters.t,2,0);
% 
%     if tlaga < 1
%         twin_cstack = [0 2];
%     else
%         twin_cstack = [tlaga-1 tlaga+1];
%     end
%     
%     [tlagc,cstackt,cstackA] = ceps_func(fader_array,Parameters.t',twin_cstack);
% 
%     if abs(tlagc - tlaga) < 0 * max(tlagc, tlaga)
%         tlag = 0.5 * (tlagc + tlaga);
%     else
%         tlag = tlaga;
%     end
% 
%     if mod(length(Parameters.t), 2)
% 
%         t           = Parameters.t(1:(length(Parameters.t)-1));
%         fader_array = fader_array(:, 1:(length(Parameters.t)-1));
% 
%     else 
% 
%         t = Parameters.t';
% 
%     end
% 
%     R_flted = filterRF_FADER(fader_array, t, tlag, r0, 0);
% 
