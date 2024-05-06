function [ ts ] = taper_ts( ts, ind1, ind2, r)
%taper_ts Put a taper on the ts with the given end points

    %make sure that the time series is horizontal vector first
    [n,~] = size(ts);
    
    if n ~= 1
        
        flip = 1;
        
        ts = ts';
       
    else
        
        flip = 0;
        
    end

    %define the initial taper
    taper = tukeywin(ind2 - ind1, r)';    
    
    %extend it to the time series length
    taper = [ zeros(1, ind1) taper zeros(1, length(ts) - ind2)];
    
    if length(taper) > length(ts)
       
        taper = taper(1:length(ts));
        
    end
    
    %apply it
    ts = ts.*taper;
        
    if flip
        
        ts = ts';
        
    end  

end

