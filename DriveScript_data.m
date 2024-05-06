%%%%no moveout correction RFs with adaptive gridding
%close all
clc
clear

warning('off', 'all')
issyn = false;

dirname    = 'TA.H55A';
%dirname    = 'Syn_Complex';

%%%%%
%surf96 has a bullshit set up requiring that parametes are loaded from a
%file with a hardwired name. So, you can't have multiple copies running. 
%so I just create multiple copies of the whole package (yes this is
%stupid). 
Parameters = make_parameters(dirname);

if exist([ dirname '_run' ]) == 7

    eval([' !rm -r ' dirname Parameters.tag '_run' ]);
    
end

mkdir([ dirname Parameters.tag '_run' ])
cd(['./' dirname Parameters.tag '_run/'])

if ~issyn

    %first step is to get the RF from EARS, so I have an event list. Yes
    %this is clunky. 

    if ~exist([ '../' dirname ])

        nm = split(dirname, '.');
        url_for_events = [ 'http://ears.iris.washington.edu/receiverFunction.zip?netCode=' nm{1} '&stacode=' nm{2} '&minPercentMatch=80&gaussian=2.5' ];
        websave('data.zip', 'http://ears.iris.washington.edu/receiverFunction.zip?netCode=TA&stacode=H55A&minPercentMatch=80&gaussian=2.5')
        !unzip data.zip
        movefile('./Ears/gauss_2.5/TA.H55A/', [ '../' dirname ])
        rmdir('Ears', 's')
        delete('data.zip')

    end

    eval([ '!cp -r ../' dirname ' ./data/'])

end

if ismac
 
    %on laptop
    eval('!cp -r ../bin_mac/ ./bin_mac/')

    path2BIN = './bin_mac/'; % path to surf96 binary
    !chmod ++x ./bin_mac/*
    addpath('../functions/')

else

    eval([ '!cp -r ../bin_linux/ ./bin_linux/' ])
    
    %on moonsoon
    path2BIN = './binlinux/';
    !chmod ++x ./bin_linux/*
    addpath('../functions-linux/')

end
%%%%%%

addpath('../toolbox/')
addpath('../rdsac/')
addpath('../deconvolution_code/')
addpath('../irisFetch/')
addpath('../Scripts/')

PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
%     setenv('PATH', [PATH,':',path2BIN]);
    setenv('PATH', [path2BIN,':',PATH]);
end

delete(gcp('nocreate'))
p                     = parpool; 
Parameters.numworkers = p.NumWorkers;

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