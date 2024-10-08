function station_string = read_station_list(filename, job_id)
    
    % Open the file
    fid = fopen(filename, 'r');
    
    % Initialize a line counter for stations
    stationLine = 1;
    
    network = {};
    station = {};

    % Loop through the file line by line
    while ~feof(fid)
        % Read the line as a string
        line = fgetl(fid);
        
        % Skip any header lines (lines starting with '#')
        if startsWith(line, '#')
            continue;
        end
        
        % Split the line by the delimiter '|'
        fields = strsplit(line, '|');
        
        % Extract the network and station name (first two fields)
        if length(fields) >= 2
            network{stationLine} = fields{1};
            station{stationLine} = fields{2};
            
            % Display the network, station, and line number
            %fprintf('Line %d: Network = %s, Station = %s\n', stationLine, network, station);
            
            % Increment station line counter
            stationLine = stationLine + 1;
        end
    end
    
    station_string = [ network{job_id} '.' station{job_id} ];

    % Close the file
    fclose(fid);

end