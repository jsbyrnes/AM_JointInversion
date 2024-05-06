function [ error ] = bootstrap_rf_amplitude( data, num, z, depth_est, depth_range)
%BOOTSTRAP Get an estimate of the error using a bootstrapping algorithm.
%This code looks for a minimum, so flip the trace to find a max. 
%Joseph Byrnes, June 2013. 

    [m, n] = size(data); %m is the length of the vectors, n is the number of samples.
    
    results = zeros(num, 1);
    
    for i = 1:num
       
        indicies = randi(n, [n 1]);
        
        tmp_distribution = zeros(size(data));
        
        tmp_distribution = data(:, indicies);
        
        tmp = mean(tmp_distribution, 2);
        
        %that sets up the bootstrap...
        %now find the depth
                
        [~, ind_plus] = findnearest(z, depth_est + depth_range/2, []);% central index
            
        [~, ind_minus] = findnearest(z, depth_est - depth_range/2, []);% central index
        
        results(i) = min(tmp(ind_minus:ind_plus));
        
        results(i) = results(i) + ind_minus;
            
    end
        
    error = 2*std(results);
    %error = std(results);
    
end

