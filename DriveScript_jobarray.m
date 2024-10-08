function DriveScript_jobarray(filename, job_id)

    %%%%no moveout correction RFs with adaptive gridding
    %close all
    
    warning('off', 'all')
    issyn = false;
    
    dirname = read_station_list(filename, job_id);

    %dirname    = 'Syn_Complex';
    Parameters = make_parameters(dirname);
        
    if issyn
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Synthetic data
        %make a new model in the format of the model structure
        rng(600)    
        truemodel.vs.z          = [ 15  15.1 22  22.1 25   35  45 150 ]';
        truemodel.vs.node       = log([ 2.5 3.5  3.5 3.0  3.5  3.5 4.5  4.5 ])';
    
        truemodel.vpvs_z     = log(0.3);
        truemodel.vpvs_block = log(4);
        truemodel.vpvs_on    = true;
    
        truemodel.max_z   = Parameters.plot_z;
    
        %uses real z
        truemodel.vs.z = truemodel.vs.z/truemodel.max_z;
        %assumes non-dimensionalized z
        truemodel.t = ((-Parameters.pre:0.05:Parameters.max_time))';
        truemodel   = update_z(truemodel, Parameters);
        truemodel   = fine_model(truemodel, Parameters);
    
        %both use real z
        truemodel.vs.z    = truemodel.vs.z*truemodel.max_z;
        Disp.c_r          = dispR_surf96(Parameters.t_array, make_surf96model(truemodel, Parameters));
        truemodel.cr_pred = Disp.c_r;
        truemodel.vs.z    = truemodel.vs.z/truemodel.max_z;
    
        %uses fine model
        rawData = make_synthetics(Parameters, truemodel, 30, 0.025);
    
        %now set other parameters needed to make a receiver function
        phase  = 'P';
        dt     = .05;
                    
        truemodel.sig_cr = -5;
        truemodel.sig_cl = -5;
        
        Disp.c_r    = Disp.c_r + 0.02*randn(size(Disp.c_r));
        Disp.c_rstd = 0.02*ones(size(Disp.c_r));
        rng('shuffle')
    
    else
        
        rng("shuffle")
        
        % Status file is located in the parent directory
        status_file = 'job_status.txt';
        
        % Check if the status file exists; if not, create it
        if ~isfile(status_file)
            fid = fopen(status_file, 'w');  % Create the file
            fclose(fid);
            disp('Job status file created.');
        end
        
        % Pause for a random amount of time between 1 and 120 seconds
        pause_time = randi([1, 120], 1, 1);
        disp(['Job ', num2str(job_id), ' pausing for ', num2str(pause_time), ' seconds...']);
        pause(pause_time);
    
        % Wait until no other job is downloading
        while true
            if ~is_any_job_downloading(status_file)
                disp(['Job ', num2str(job_id), ' starting downloads...']);
                break;
            end
            disp(['Job ', num2str(job_id), ' is waiting for other jobs to finish downloading...']);
            pause(60);  % Wait for 60 seconds before checking again
        end
    
        % Mark the current job as downloading
        mark_job_as_downloading(job_id, status_file);
    
        % Execute ConfigureRun located in the parent directory
        run('ConfigureRun');
    
        % Run the data download process
        disp('-----> Getting the data...');
        [rawData, Disp, Parameters] = load_real_data(dirname, Parameters);
    
        % After the download is complete, mark the job as finished
        mark_job_as_finished(job_id, ['../'  status_file ]);
    
        truemodel = [];
    
    end
    
    Parameters.t = ((-Parameters.pre:0.05:Parameters.max_time))';
    
    %%%%%%%%
    %Run the inverse problem
    [ model, allWfs ] = iterative_inversion(Parameters, rawData, Disp, truemodel);
    
    %%%%%%%%%%%%%%
    %Bootstrapping
    if Parameters.bootstrapping
    
        model = bootstrapping(model, Parameters, rawData, Disp, truemodel);
    
    end
    
    save([ Parameters.save_name 'inversion.mat'], 'Parameters', 'allWfs', 'model', 'Disp', 'truemodel', 'rawData')
    %%%%%%%%

end

function downloading = is_any_job_downloading(status_file)
    % Check if the status file indicates an ongoing download or if it's empty (indicating nothing has started).
    
    downloading = false;  % Assume no jobs are downloading initially
    
    if exist(status_file, 'file')
        fid = fopen(status_file, 'r');
        lastLine = '';  % Variable to store the last line
        isEmpty = true;  % Flag to check if the file is empty
        
        % Read through the file to find the last line
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line)
                lastLine = line;  % Store the current line as the last one
                isEmpty = false;  % File has content
            end
        end
        
        % Check if the file was empty
        if isEmpty
            downloading = false;  % If file is empty, assume nothing has started
        else
            % Check the last line for the specific status
            if contains(lastLine, 'finished downloading')
                downloading = false;
            elseif contains(lastLine, 'downloading')
                downloading = true;
            end
        end
        
        fclose(fid);
    else
        % If the file does not exist, assume nothing has started
        downloading = false;
    end
end

function mark_job_as_downloading(job_id, status_file)
    % Append to the status file that this job is downloading
    fid = fopen(status_file, 'a');
    fprintf(fid, 'Job_id %d downloading\n', job_id);
    fclose(fid);
end

function mark_job_as_finished(job_id, status_file)
    % Append to the status file that this job has finished downloading
    fid = fopen(status_file, 'a');
    fprintf(fid, 'Job_id %d finished downloading\n', job_id);
    fclose(fid);
end
