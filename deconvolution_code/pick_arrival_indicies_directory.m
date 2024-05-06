%for a given events station pair you need to know where to 'clip' the trace
%for RFs analysis. This is done by ginput in this script.

%To help guide you, the S (and SKS if it is there) arrival times are
%plotted as verical lines.
%Then you select the edges of a taper to be applied during deconvolution to
%the R trace.

close all

directory = '../CI_yr1_PS/';

phase = 'P';

r = 0.3;

if strcmp('SV', phase)
    
    str = 'R trace';
    
elseif strcmp('SH', phase)
    
    str = 'T trace';
    
elseif strcmp('P', phase)
    
    str = 'R';
    
end

files = dir(directory);
files = {files.name};
files = files(3:end);

for j = 1:length(files)
    
    disp([ 'On file ' files{j} ', number ' num2str(j) ' of ' num2str(length(files)) ]);
    load([ directory files{j} ]);
    
    for i = 1:length(clTrace)
        
        disp([ 'On trace ' num2str(i) ' of ' num2str(length(clTrace)) ]);
        
        if (isempty(clTrace(i).dt))
            
            continue
            
        end
        
        %the predicted arrival time is at time 0
        t = ((1:length(clTrace(i).Z_trace))*clTrace(i).dt) - ((length(clTrace(i).Z_trace) - 1)/2)*clTrace(i).dt;
        
        if strcmp('SV', phase) || strcmp('P', phase)
            
            plot_trace = clTrace(i).R_trace;
            
        elseif strcmp('SH', phase)
            
            plot_trace = clTrace(i).T_trace;
            
        end
        
        if isempty(clTrace(i).Z_trace)
            
            continue
            
        end
        
        if isempty(plot_trace)
            
            continue
            
        end
        
        figure(1)
        
        good = 'n';
        
        while(good(1) ~= 'y')
            
            clf
            
            hold on
            plot(t, plot_trace);
            plot(t, clTrace(i).Z_trace);
            legend(str,'Z trace');
            xlabel('Time (seconds)');
            ylabel('Amplitude (in counts)');
            legend(str, 'Z trace');
            
            [ t_picked, ~ ] = ginput(2);
            
            t_picked = sort(t_picked);
            
            [~, ind1] = min(abs(t - t_picked(1)));
            [~, ind2] = min(abs(t - t_picked(2)));
            
            if strcmp('P', phase)
            
                trace_new = taper_ts(clTrace(i).Z_trace, ind1, ind2, r);
            
                clf

                hold on
                plot(t, plot_trace);
                plot(t, trace_new);
                xlabel('Time (seconds)');
                ylabel('Amplitude (in counts)');
                legend(str, 'Z trace');
                
            else
               
                trace_new = taper_ts(plot_trace, ind1, ind2, r);
            
                clf

                hold on
                plot(t, trace_new);
                plot(t, clTrace(i).Z_trace);
                xlabel('Time (seconds)');
                ylabel('Amplitude (in counts)');
                legend(str, 'Z trace');
                
            end
            
            y_limits = ylim;
            
            line([ t_picked(1) t_picked(1) ], [-1e10 1e10], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', [0 0 0]);
            line([ t_picked(2) t_picked(2) ], [-1e10 1e10], 'LineStyle', '--', 'LineWidth', 0.5, 'Color', [0 0 0]);
            
            ylim(y_limits);
            
            good = input('Good?', 's');
            
            if length(good) ~= 1
                
                good = 'n';
                
            end
            
        end
        
        clTrace(i).taper_indicies = [ ind1 ind2 ];
        
    end
    
    save([ directory files{j} ], 'clTrace');
    
end
