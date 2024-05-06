function [ error ] = bootstrap_mean( data, num )
%BOOTSTRAP Get an estimate of the error using a bootstrapping algorithm.
%data in columns, row are seperate measurements. num is the number of
%bootstrap samples to make. Returns the bootstrap error at each point along
%the vector (m by 1 returned)
%Joseph Byrnes, June 2013. 

    [m, n] = size(data); %m is the length of the vectors, n is the number of samples.
    
    results = zeros(m, num);
    
    for i = 1:num
       
        indicies = randi(n, [n 1]);
        
        tmp_distribution = zeros(size(data));
        
        tmp_distribution = data(:, indicies);
        
        tmp = mean(tmp_distribution, 2);
        
        results(:, i) = tmp;
            
    end
        
    error = std(results, 0, 2);

end

