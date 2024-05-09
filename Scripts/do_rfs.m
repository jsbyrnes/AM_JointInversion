function synWfs = do_rfs(model, Parameters, allWfs, rotate)

    if nargin==3

        rotate = true;
        
    end

    synWfs = allWfs;

    [~, G] = size(allWfs);

    %sue me
    model.vs   = model.vs.model;

    model.vpvs = model.vp./model.vs;

%     if isfield(model, 'vpvs')
% 
%         model.vpvs = model.vpvs.model;
% 
%     else
%        
%         model.vpvs = 1.76*ones(size(model.vs));
% 
%     end

    [P_comp, SV_comp, SH_comp] = anirec("P", 0.05, Parameters.move_out.ref_p, 0, model, 0);
    [~, ind] = max(P_comp);

    indn    = round(Parameters.pre/0.05) + 1;
    P_comp  = circshift(P_comp, indn - ind);
    SV_comp = circshift(SV_comp, indn - ind);
    SH_comp = circshift(SH_comp, indn - ind);

    f = @(tau) centering(P_comp, 1/0.05, tau, indn);

    tau = fminbnd(f, -0.05, 0.05);

    P_comp  = delay_continuous(P_comp, 1/0.05, tau);
    SV_comp = delay_continuous(SV_comp, 1/0.05, tau);

    if rotate

        [ P_comp, SV_comp ] = ZR2PSV( P_comp, -1*SV_comp, model.vp(1), model.vp(1)/model.vpvs(1), Parameters.move_out.ref_p);
        P_comp  = -1*P_comp;
        SV_comp = -1*SV_comp;

    end

    P_comp  = [ P_comp; zeros(size(P_comp)) ];
    SV_comp = [ SV_comp; zeros(size(SV_comp)) ];
    SH_comp = [ SH_comp; zeros(size(SH_comp)) ];

    %[E,V] = dpss(800, 4, 7);

    for j = 1:G

        window = max([ 3200 - 800*(j-1), 800]);
        
        %window = 400;
        [E,V] = dpss(window/8, 2.5, 3);

        if Parameters.deconvolve 

            if strcmp(Parameters.rf_style, 'multi-taper')

                [~, rf, ~, ~] = multitaper2rf_3component(P_comp(1:window), SV_comp(1:window), ...
                    SH_comp(1:window), 1/Parameters.sample_rate, Parameters.pre, E,V, 'P', ...
                    [ 0.05 Parameters.g_array(j) ], [ 1 length(P_comp)]);

            elseif strcmp(Parameters.rf_style, 'time-domain')

                rf = make_RF(P_comp, SV_comp, 0.05, Parameters.g_array(j), Parameters.pre, Parameters.filter_style);

            end

        else

            if strcmp(Parameters.filter_style, 'gaussian')

                gsig = 1/(pi*Parameters.g_array(j));
                gauss_sig = round(gsig/0.05);

                x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
                Gauss_win = exp(-x.^2/(2*gauss_sig^2));
                P = conv(P_comp,Gauss_win,'same');
                D = conv(SV_comp,Gauss_win,'same');

                rf = D/max(abs(P));

            elseif strcmp(Parameters.filter_style, 'butter')

                P = bandpassfilt_rfs(P_comp, 1/Parameters.sample_rate, Parameters.g_array(j), 0.05);
                D = bandpassfilt_rfs(SV_comp, 1/Parameters.sample_rate, Parameters.g_array(j), 0.05);                

                rf = D/max(abs(P));

            end

        end
        
        t  = ((0:length(rf)-1))*0.05 - Parameters.pre;

        rf = interp1(t, rf, model.t, 'linear');

        if isrow(rf)

            rf = rf';

        end

        synWfs(1,j).rfr = rf.*tukeywin(length(rf), 0.1);

    end

end


%     parfor k = 1:length(Parameters.p_array)
% 
%         slow = Parameters.p_array(k);
% 
%         [P_comp, SV_comp, SH_comp] = anirec("P", 0.05, slow, 0, model, 0);
%         [~, ind] = max(P_comp);
% 
%         indn    = round(Parameters.pre/0.05) + 1;
%         P_comp  = circshift(P_comp, indn - ind);
%         SV_comp = circshift(SV_comp, indn - ind);
%         SH_comp = circshift(SH_comp, indn - ind);
% 
%         f = @(tau) centering(P_comp, 1/0.05, tau, indn);
% 
%         tau = fminbnd(f, -0.05, 0.05);
% 
%         P_comp  = delay_continuous(P_comp, 1/0.05, tau);
%         SV_comp = delay_continuous(SV_comp, 1/0.05, tau);
% 
%         [ P_comp, SV_comp ] = ZR2PSV( P_comp, -1*SV_comp, model.vp(1), model.vp(1)/model.vpvs(1), slow);
%         P_comp  = -1*P_comp;
%         SV_comp = -1*SV_comp;
% 
%         for j = 1:G
%     
%             if Parameters.deconvolve 
% 
%                 if strcmp(Parameters.rf_style, 'multi-taper')
% 
%                     [~, rf, ~] = multitaper2rf_3component(P_comp(1:2000), SV_comp(1:2000), ...
%                         SH_comp(1:2000), 1/Parameters.sample_rate, Parameters.pre, 400, 4, 7, 'P', ...
%                         [ 0.05 Parameters.g_array(j) ], [ 1 length(P_comp)]);
% 
%                 elseif strcmp(Parameters.rf_style, 'time-domain')
% 
%                     rf = make_RF(P_comp, SV_comp, 0.05, Parameters.g_array(j), Parameters.pre, Parameters.filter_style);
% 
%                 end
% 
%             else
% 
%                 if strcmp(Parameters.filter_style, 'gaussian')
% 
%                     gsig = 1/(pi*Parameters.g_array(j));
%                     gauss_sig = round(gsig/0.05);
% 
%                     x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
%                     Gauss_win = exp(-x.^2/(2*gauss_sig^2));
%                     P = conv(P_comp,Gauss_win,'same');
%                     D = conv(SV_comp,Gauss_win,'same');
% 
%                     rf = D/max(abs(P));
% 
%                 elseif strcmp(Parameters.filter_style, 'butter')
% 
%                     P = bandpassfilt_rfs(P_comp, 1/Parameters.sample_rate, Parameters.g_array(j), 0.05);
%                     D = bandpassfilt_rfs(SV_comp, 1/Parameters.sample_rate, Parameters.g_array(j), 0.05);                
% 
%                     rf = D/max(abs(P));
% 
%                 end
% 
%             end
%             
%             t  = ((0:length(rf)-1))*0.05 - Parameters.pre;
% 
%             rf = interp1(t, rf, model.t, 'linear');
% 
%             if isrow(rf)
% 
%                 rf = rf';
% 
%             end
% 
%             synWfs(k,j).rfr = rf.*tukeywin(length(rf), 0.2);
%     
%         end
% 
%     end

