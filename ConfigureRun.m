%This script contains common calls that set up all the runs. 
%Do not run this by itself or modify unless you have a good
%reason. 

%%%%%

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
        evalc('!unzip data.zip')
        movefile([ './Ears/gauss_2.5/' dirname '/' ], [ '../' dirname ])
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
    path2BIN = './bin_linux/';
    !chmod ++x ./bin_linux/*
    addpath('../functions-linux/')

end

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

%%%%%%
