function tmp = plot_model_new(number, model, Parameters, allWfs, Disp, truemodel)

    if ~Parameters.plot
        tmp = 0;
        return
    end

    panels = 2 + length(Parameters.fields_vec);

    model_fine = fine_model(model, Parameters);

    figure(number)
    clf

    for k = 1:length(Parameters.fields_vec)

        field = Parameters.fields_vec{k};

        subplot(1,panels,k)
        if ~isempty(truemodel)    
    
            truemodel = fine_model(truemodel, Parameters);
            
            if truemodel.z(1) > 0

                truemodel.z = [ 0; truemodel.z ];
                truemodel.(field).model = [ truemodel.(field).model(1); truemodel.(field).model ];

            end

            stairs(truemodel.(field).model, truemodel.z, 'r', 'LineWidth', 1)
            hold on
            
        end    
        
        if model_fine.z(1) == 0
    
            model_fine.z =  model_fine.z(2:end);
            model_fine.(field).model = model_fine.(field).model(1:length(model_fine.z));
    
            model_fine.z = [ 0; model_fine.z];
            model_fine.(field).model = [ model_fine.(field).model(1); model_fine.(field).model];

        end

        stairs(model_fine.(field).model, model_fine.z, 'k', 'LineWidth', 2)%[ 68, 119, 170 ]/256
        hold on%ok if used twice, but usually bad if used before any plotting
        plot(exp(model.(field).node), model.(field).z*model.max_z, 'ko')
        set(gca, 'YDir', 'reverse')
        xlabel(field)
        ylabel('Depth, km')
        ylim([ 0 Parameters.plot_z ])
        %xl = xlim;
    
        if strcmp(field, 'vs')

            xticks([ 1 2 3:0.5:5 ])

        end

        if model.vpvs_on && k == 1

            title([ 'VpVs of ' num2str(exp(model.vpvs_block)) ' to ' num2str(exp(model.vpvs_z)) ' km'])

        end

        %if diff(xl) < 1
    
        %    xlim(mean(xl) + [ -1.5 1.5 ]);
    
        %end
    
        grid on

    end

    subplot(1,panels,length(Parameters.fields_vec) + 1)    
    plot(Parameters.t_array, model.cr_pred, 'r', 'LineWidth', 2)
    hold on
    errorbar(Parameters.t_array, Disp.c_r, Disp.c_rstd, 'ko')
    set(gca, 'XScale', 'log')    
    xlabel('Period, s')
    ylabel('Phase velocity, km/s')
    grid on

    [p,g] = size(allWfs);

    subplot(1,panels,length(Parameters.fields_vec) + 2)
    hold on
    for k = 1:g

        plot(model.t, 0.5*allWfs(p, k).rfr/max(abs(allWfs(p, k).rfr)) + k, 'r', 'LineWidth', 1)
        %plot(model.t, (0.5*(allWfs(p, k).rfr + allWfs(p, k).rfr_std))/max(abs(allWfs(p, k).rfr)) + k, 'r--')
        %plot(model.t, (0.5*(allWfs(p, k).rfr - allWfs(p, k).rfr_std))/max(abs(allWfs(p, k).rfr)) + k, 'r--')
        plot(model.t, 0.5*model.synWfs(p, k).rfr/max(abs(allWfs(p, k).rfr)) + k, 'k', 'LineWidth', 1)

        yl{k} = num2str(Parameters.g_array(k));

    end

    yticks(1:length(Parameters.g_array(1:g)));
    yticklabels(yl)
    ylabel('Filter Frequency, Hz')
    xlabel('Time, s')

    drawnow
    
    print(number, [ './' Parameters.save_name 'Model.pdf' ], '-dpdf', '-bestfit')

    tmp = 1;

end

%     subplot(1,panels,2)
%     if ~isempty(truemodel)  
%         
%         stairs([ (truemodel.vpvs.model(1)); (truemodel.vpvs.model) ], [ 0; truemodel.z ], 'r', 'LineWidth', 2')
%         hold on
% 
%     end    
%     
%     stairs([ (model_fine.vpvs.model(1)); (model_fine.vpvs.model) ], [ 0; model_fine.z ], 'k', 'LineWidth', 2')
%     set(gca, 'YDir', 'reverse')
%     xlabel('Vp/Vs')
%     ylabel('Depth, km')
%     ylim([ 0 Parameters.plot_z ])
%     grid on
% 
%     xl = xlim;
% 
%     if diff(xl) < 0.2
% 
%         xlim(mean(xl) + [ -0.1 0.1 ]);
% 
%     end
