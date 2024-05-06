function model = bootstrapping(model, Parameters, rawData, Disp, truemodel)

    p0 = Parameters.p_array;

    for k = 1:Parameters.bts_samples

        %gets modified by bootstrapping sometimes
        Parameters.p_array = p0;

        disp([ '-----> Bootstrap sample #' num2str(k) '...' ])

        if ~isempty(Parameters.baz)
    
            if Parameters.baz(2) > Parameters.baz(1)
    
                ind_set = [rawData(:).baz]>Parameters.baz(1) & [rawData(:).baz] < Parameters.baz(2);
    
            else
    
                ind_set = [rawData(:).baz]>Parameters.baz(2) & [rawData(:).baz] < Parameters.baz(1);
    
            end
    
            ind     = randi(length(ind_set), length(ind_set), 1);

        else
    
            ind     = randi(length(rawData), length(rawData), 1);
    
        end

    	allWfs  = fill_RFs(model(1).vs.model(1), model(1).vs.model(1)*exp(model(1).vpvs_block), Parameters.t, rawData(ind), Parameters, 500);

        %some get skipped!
        kill = zeros(length(allWfs(:, 1)),1);
        
        for ix = 1:length(Parameters.p_array)
        
            if isempty(allWfs(ix, 1).rfr)
        
                kill(ix) = 1;
        
            end
        
        end
        
        if any(kill)
        
            allWfs(kill==1, :)         = [];
            Parameters.p_array(kill==1) = [];
        
        end

        Dispn     = Disp;
        Dispn.c_r = normrnd(model(1).cr_pred, Disp.c_rstd);

        mn = evaluate_reflectivity(model(1), allWfs, Dispn, Parameters, false);

        %start at the optimum
        %model(1+k) = inversion(model(1+k), Parameters, allWfs, Dispn, truemodel);

        %scramble the grid a little
        mn = expand_grid(mn, Parameters, Disp, allWfs, 0.5);
        model(1+k) = inversion(mn, Parameters, allWfs, Dispn, truemodel);

    end

    %and get model stats

    z = 0:0.1:Parameters.plot_z;

    for k = 1:length(z)
    
        for j = 2:length(model)
    
            vs_bts(j-1) = interp1(model(j).z, model(j).vs.model, z(k), Parameters.interp, model(j).vs.model(1));

        end
    
        vs_m(k)   = mean(vs_bts);
        vs_std(k) = std(vs_bts);
    
    end
    
    figure(1)
    if ~isempty(truemodel)    
    
        truemodel = fine_model(truemodel, Parameters);
        
        stairs([ truemodel.vs.model(1); truemodel.vs.model ], [ 0; truemodel.z ], 'r', 'LineWidth', 2')
        hold on
        
    end    
    
    figure(1)
    stairs([ vs_m(1); vs_m' ], [ 0; z' ], 'k', 'LineWidth', 2)
    hold on
    stairs([ vs_m(1) + vs_std(1); (vs_m + vs_std)' ], [ 0; z' ], 'k--', 'LineWidth', 1)
    stairs([ vs_m(1) - vs_std(1); (vs_m - vs_std)' ], [ 0; z' ], 'k--', 'LineWidth', 1)
    stairs([ vs_m(1) + 2*vs_std(1); (vs_m + 2*vs_std)' ], [ 0; z' ], 'k--', 'LineWidth', 0.25)
    stairs([ vs_m(1) - 2*vs_std(1); (vs_m - 2*vs_std)' ], [ 0; z' ], 'k--', 'LineWidth', 0.25)
    set(gca, 'YDir', 'reverse')
    xlabel('Vs, km/s')
    ylabel('Depth, km')
    ylim([ 0 Parameters.plot_z ])
    xl = xlim;
    grid on
    title('Bootstrapped Result')
    print(1, 'BootstrappedResult.pdf', '-dpdf')

end