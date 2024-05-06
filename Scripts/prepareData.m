function allWfs = prepareData(vs_surface, t_interp, rawData, Parameters)

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

    for k = 1:length(rawData)
        
        [ P_comp, SV_comp ] = ZR2PSV(rawData(k).Z, -1*rawData(k).R, ...
            vs_surface*2, vs_surface, rawData(k).p);%surface vs won't matter
        P_comp  = -1*P_comp;
        SV_comp = -1*SV_comp;
        SH_comp = zeros(size(SV_comp));%need to flip?
    
        for j = 1:G
                
            filterlims = [ 0.01 Parameters.g_array(j) ];
            
            P = bandpassfilt_rfs(P_comp, 1/Parameters.sample_rate, filterlims(2), filterlims(1));
            D = bandpassfilt_rfs(SV_comp, 1/Parameters.sample_rate, filterlims(2), filterlims(1));                

            t   = ((0:length(P)-1))*0.05 - Parameters.pre;
            
            P = interp1(t, P, t_interp, 'linear');
            D = interp1(t, D, t_interp, 'linear');
            
            if isrow(P)
            
                P = P';
                D = D';
            
            end
                            
            allWfs(k,j).Z = P;
            allWfs(k,j).R = D;

            allWfs(k,j).p = rawData(k).p;
            allWfs(k,j).g = Parameters.g_array(j);
                
        end

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
