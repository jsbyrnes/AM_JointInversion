%%%%no moveout correction RFs with adaptive gridding
%close all
clc
clear

warning('off', 'all')
issyn = true;

%dirname    = 'TA.H55A';
dirname    = 'Syn_Complex';
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
    truemodel.vpvs_block = log(1.76);
    truemodel.vpvs_on    = false;

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
    rawData = make_synthetics(Parameters, truemodel, 30, 0.0001);

    %now set other parameters needed to make a receiver function
    phase  = 'P';
    dt     = .05;
                
    truemodel.sig_cr = -5;
    truemodel.sig_cl = -5;
    
    Disp.c_r    = Disp.c_r + 0.025*randn(size(Disp.c_r));
    Disp.c_rstd = 0.02*ones(size(Disp.c_r));
    rng('shuffle')

else

    disp('-----> Getting the data...')
    [rawData, Disp, Parameters] = load_real_data(dirname, Parameters);

    truemodel = [];

end

Parameters.t = ((-Parameters.pre:0.05:Parameters.max_time))';

%%%%%%%%
%Run the inverse problem
[ model, allWfs ] = iterative_inversion(Parameters, rawData, Disp, truemodel);
%[ model, allWfs ] = node_inversion(Parameters, rawData, Disp, truemodel, []);

%%%%%%%%%%%%%%
%Bootstrapping
if Parameters.bootstrapping

    model = bootstrapping(model, Parameters, rawData, Disp, truemodel);

end

save([ Parameters.save_name 'inversion.mat'], 'Parameters', 'allWfs', 'model', 'Disp', 'truemodel', 'rawData')
%%%%%%%%