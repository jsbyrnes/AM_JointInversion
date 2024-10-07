function DriveScript_jobarray(filename, job_id)

    %%%%no moveout correction RFs with adaptive gridding
    %close all
    
    % Define file to track the job statuses
    status_file = 'job_status.txt'; %This is used internally to prevent the job arrays from downloading at the same time. 

    warning('off', 'all')
    issyn = false;
    
    dirname = read_station_list(filename, job_id);

    %dirname    = 'Syn_Complex';
    Parameters = make_parameters(dirname);
    
    ConfigureRun
    
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
    
        % If job_id is 1, create the file and leave it blank
        if job_id == 1
            fid = fopen(status_file, 'w');
            fclose(fid);  % Blank file is created
        else
            prev_job_id = job_id - 1;
            disp(['Job ', num2str(job_id), ' is waiting for job ', num2str(prev_job_id), ' to finish downloading...']);
    
            % Keep checking the status file until the previous job is done
            while true
                if is_job_finished(prev_job_id, status_file)
                    disp(['Job ', num2str(prev_job_id), ' finished. Job ', num2str(job_id), ' starting downloads...']);
                    break;
                end
                pause(60);  % Wait for 60 seconds before checking again
            end
        end
        
        % Run the data download process
        disp('-----> Getting the data...');
        [rawData, Disp, Parameters] = load_real_data(dirname, Parameters);
    
        % After the download is complete, write to the status file
        fid = fopen(status_file, 'a');
        fprintf(fid, 'Job_id %d finished downloading\n', job_id);
        fclose(fid);
    
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

function finished = is_job_finished(job_id, status_file)
    % Check if the status file contains the "finished" message for the given job_id
    finished = false;
    if exist(status_file, 'file')
        fid = fopen(status_file, 'r');
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, sprintf('Job_id %d downloading... finished', job_id))
                finished = true;
                break;
            end
        end
        fclose(fid);
    end
end
