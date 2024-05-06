function rawData = make_synthetics(Parameters, truemodel, nsyn, noise)
    
    Parameters.p_array = rand(nsyn, 1)*0.04 + 0.04;
    rawData            = do_syns(truemodel, Parameters);

    t = (0.05*(1:length(rawData(1).Z)) - Parameters.pre)';

    for k = 1:length(Parameters.p_array)

        if Parameters.deconvolve

            %how many sources?
    
            ns = randi(1000) + 1000;
    
            source = zeros(size(rawData(k).Z));
    
            for si = 1:ns
    
                order = ceil(exprnd(10));
                shift = round((exprnd(8))/0.05);
                %sigma = min([ max([ exp(normrnd(0,2)) 0.05 ]) 5 ]);
                sigma = rand()*(5 - 0.025) + 0.025;
                amp   = normrnd(0, 2);
    
                s           = ((t/sigma).^2.*exp(-(t/sigma).^2)).*polyval(HermitePoly(order), t/sigma);
                s(t<0)      = 0;
                s(isnan(s)) = 0;
                s           = circshift(s, shift);
                source      = source + amp*s/sum(s.^2);
               
            end

            %undo the pre, get's double applied
            source = circshift(source, -1*Parameters.pre/0.05);

            rawData(k).Z = ifft(fft(rawData(k).Z).*fft(source), 'symmetric');
            rawData(k).R = ifft(fft(rawData(k).R).*fft(source), 'symmetric');

        end

        %tpicking = (0:(round(length(rawData(k).Z))-1))/Parameters.sample_rate - Parameters.pre;

%         [pick, snr] = sta_lta(rawData(k).Z, -Parameters.pre, 50, ...
%             sta, lta, tpicking);
%     
%         [~, ind] = min(abs(tpicking - pick));
    
        ind  = Parameters.sample_rate*Parameters.pre;
        ind2 = abs(Parameters.datawin(1))*Parameters.sample_rate;
        %only keep 100 s
        try
    
            rawData(k).Z = [ zeros(ind2, 1); rawData(k).Z ];
            rawData(k).R = [ zeros(ind2, 1); rawData(k).R ];

            rawData(k).Z = rawData(k).Z((ind + ind2 - abs(Parameters.datawin(1))*Parameters.sample_rate):(ind + 200*Parameters.sample_rate));
            rawData(k).R = rawData(k).R((ind + ind2 - abs(Parameters.datawin(1))*Parameters.sample_rate):(ind + 200*Parameters.sample_rate));
    
        catch
    
            disp(' -> Trace length incorrect. Pick probably too early')
            continue%very rare
    
        end


        rawData(k).R = rawData(k).R/max(abs(rawData(k).Z));
        rawData(k).Z = rawData(k).Z/max(abs(rawData(k).Z));
        
        rawData(k).Z = rawData(k).Z + noise*randn(size(rawData(k).Z));
        rawData(k).R = rawData(k).R + noise*randn(size(rawData(k).Z));
    
        taper = (abs(Parameters.datawin(1))/(0.5*(0.05*length(rawData(k).Z))));
        rawData(k).Z = rawData(k).Z.*tukeywin(length(rawData(k).Z), taper);
        rawData(k).R = rawData(k).R.*tukeywin(length(rawData(k).R), taper);

    end

end