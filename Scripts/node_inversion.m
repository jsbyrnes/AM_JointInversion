function model = node_inversion(Parameters, allWfs, Disp, truemodel, model0)

    if nargin == 4 || isempty(model0)

        %%%%%%%%%    
        disp('-----> Building a starting model')        
        model.vs.z      = cumsum(1:1:10)';
        model.vs.z      = unique([ model.vs.z; (model.vs.z(end):10:Parameters.plot_z)' ]);
        model.vs.z(model.vs.z>Parameters.plot_z) = [];
        model.vs.z(end) = Parameters.plot_z;
        model.vs.z      = model.vs.z/Parameters.plot_z;
            
        model.vs.n      = length(model.vs.z);
        %model.vs.node   = log(3, 4, length(model.vs.z))';
        model.vs.node   = log(linspace(min(Disp.c_r), max(Disp.c_r), length(model.vs.z)))';%log((5 - 2.^linspace(2,0,length(model.vs.z)))');
    
        if any(contains(Parameters.fields_vec, 'vpvs'))
    
            model.vpvs.z      = cumsum(0:1:20)';
            model.vpvs.z      = unique([ model.vpvs.z; (model.vpvs.z(end):20:Parameters.plot_z)' ]);
            model.vpvs.z(model.vpvs.z>Parameters.plot_z) = [];
            model.vpvs.z(end) = Parameters.plot_z;
            model.vpvs.z      = model.vpvs.z/Parameters.plot_z;
            model.vpvs.node   = 1.8*ones(size(model.vpvs.z));
            model.vpvs.n      = length(model.vpvs.z);
    
        end
    
        model.sig_rf     = log(0.5)*ones(length(Parameters.g_array), 1);
        model.max_z      = Parameters.plot_z;
        model            = update_z(model, Parameters);

        %always defined, may not be used interrnally
        model.vpvs_z     = log(0.25);%log
        model.vpvs_block = log(1.76);%log
        model.vpvs_on         = false;%track here, not with Parameters

        %model.vs.z = model.vs.z*model.max_z;%?
        %model.z = (0:Parameters.dz_lim(1):Parameters.plot_z)';
        model      = fine_model(model, Parameters);
        %model.vs.z = model.vs.z/model.max_z;%?
    
        %allWfs         = fill_RFs(model.vs.model(1), model.vs.model(1)*exp(model.vpvs_block), ...
        %    Parameters.t, rawData, Parameters, 1);

        model.t = ((1:length(allWfs(1).rfr))*0.05 - Parameters.pre)';

        Parameters.anirec          = true;
        model                      = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
        
        plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
        %%%%%%%%
       
        %%%%%%%%
        disp('-----> Initializing to the surface waves')
        Parameters.anirec     = false;
        model                 = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
        model                 = inversion(model, Parameters, allWfs, Disp, truemodel);
        Parameters.anirec     = true;
        model                 = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
        plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
        
        if ~Parameters.move_out.perform_moveout
        
            plot_data(101, Parameters, model.synWfs, model.Disp);
    
        end

        %test for an ultra layer
        %[ac, l] = xcorr(allWfs(end).rfr, allWfs(end).rfr, 'coeff');
        %x = findpeaks(ac(l>0));

        t    = Parameters.t(Parameters.t>0 & Parameters.t<4);
        rfr  = allWfs(end).rfr(Parameters.t>0 & Parameters.t<=4);
        x    = (abs(hilbert(rfr)));
        x    = x(1:end-1);%last sample can be weird
        func = @(y) sum((y(2)*exp(-t/y(1)) + y(3) - x).^2);
        y    = particleswarm(func, 3, [ 0.01 0 0 ], [ 4 1 1]);
 
        if y(1) > 1 && y(2) > 0.1
            
            disp('Adding a fine sediment layer')
            %model = add_sed_layer(model, allWfs, Disp, Parameters, truemodel);

            %add it in to the model, a random guess
            model.vs.z       = [ ([ 0.05; 0.1; 0.2  ]/Parameters.plot_z);       model.vs.z ];
            model.vs.node    = [ log([ 0.05; 0.2; 0.5 ]*exp(model.vs.node(1))); model.vs.node ];
            model.vpvs_on    = true;
            model.vpvs_block = log(3);

            model  = update_z(model, Parameters);
            
        end

    else

        model = model0;
        
        model         = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
        plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);

    end
    %%%%%%%%
    model.t       = ((1:length(allWfs(1).rfr))*0.05 - Parameters.pre)';
    model.sig_rf  = log(0.05)*ones(length(Parameters.g_array), 1);%not used anymore
    model         = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
    
    plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);

    if ~Parameters.move_out.perform_moveout

        plot_data(101, Parameters, model.synWfs, model.Disp);

    end
        
    disp('First fit with RFs')
    model = inversion(model, Parameters, allWfs, Disp, truemodel);
    model = reduce_grid(model, Parameters, Disp, allWfs, truemodel);
    plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
    
    if ~Parameters.move_out.perform_moveout
    
        plot_data(101, Parameters, model.synWfs, model.Disp);

    end

    keepgoing = true;
    while keepgoing

        model = update_z(model, Parameters);
        model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);

        IC0 = model.IC;
        
        disp('Expanding the grid')
        mn = expand_grid(model, Parameters, Disp, allWfs);

        model = mn;

        model = inversion(model, Parameters, allWfs, Disp, truemodel);
        
        %try to change vpvs on the expanded model 
        if ~model.vpvs_on
    
            disp('Testing Vp/Vs variations')
            model.vpvs_on = true;
            model2         = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
            model2         = inversion(model2, Parameters, allWfs, Disp, truemodel);
    
            if model2.IC > model.IC
                
                disp('Keeping Vp/Vs variations')
                model = model2;
    
            else
    
                disp('Rejecting Vp/Vs variations')
                model.vpvs_on = false;
                
            end
    
        else
    
            disp('Testing no Vp/Vs variations')
            model.vpvs_on = false;
            model2         = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
            model2         = inversion(model2, Parameters, allWfs, Disp, truemodel);
    
            if model2.IC > model.IC
                
                disp('Removing Vp/Vs variations')
                model = model2;
        	    model.vpvs_z     = log(0.25);%log
	            model.vpvs_block = log(1.76);%log
    
            else
    
                disp('Keeping Vp/Vs variations')
                model.vpvs_on = true;
                
            end
    
        end        
        
        model = reduce_grid(model, Parameters, Disp, allWfs, truemodel);
        model = inversion(model, Parameters, allWfs, Disp, truemodel);
        dl = model.IC - IC0;

        disp([ ' --> Refined for an ' Parameters.IC ' of ' num2str(model.IC), ...
            '; improvement after z update of ' num2str(dl) ])

        plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);

        if ~Parameters.move_out.perform_moveout

            plot_data(101, Parameters, model.synWfs, model.Disp);

        end

        if dl < 1 %move on with microimprovements

            keepgoing = false;

        end
        
    end

end

%     keepgoing = true;
%     while keepgoing
% 
%         model = update_z(model, Parameters);
%         model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
% 
%         IC0 = model.IC;
% 
%         model = update_grid(model, Parameters, Disp, allWfs, truemodel);
%         model = inversion(model, Parameters, allWfs, Disp, truemodel);
% 
%         dl = model.IC - IC0;
% 
%         disp([ ' --> Refined for an ' Parameters.IC ' of ' num2str(model.IC), ...
%             '; improvement after z update of ' num2str(dl) ])
% 
%         if dl < 0
% 
%             disp('     Fit can degrade when updating the z vector')
% 
%         end
%         
%         plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
%         plot_data(101, Parameters, model.synWfs, model.Disp);
% 
%         if dl < 1 && dl >= 0 %don't quit if its a z-update problem degrading the fit
% 
%             keepgoing = false;
% 
%         end
%         
%     end


%     g0 = Parameters.g_array;
%         
%     model.sig_rf     = log(0.05);
% 
%     for idx = 1:length(g0)
%     
%         Parameters.g_array = g0(1:idx);
%         disp([ '-----> Adding in ' num2str(Parameters.g_array(idx)) ' Hz body waves' ])
%     
%         allWfsn = allWfs(:, 1:idx);
%         model = update_z(model, Parameters);
%     
%         model   = evaluate_reflectivity(model, allWfsn, Disp, Parameters, false);
%         model   = inversion(model, Parameters, allWfsn, Disp, truemodel);
%         plot_model_new(5678, model, Parameters, allWfsn, Disp, truemodel);
%         plot_data(101, Parameters, model.synWfs, model.Disp);
%     
%         keepgoing = true;
%         while keepgoing
%     
%             model = update_z(model, Parameters);
%             model = evaluate_reflectivity(model, allWfsn, Disp, Parameters, false);
%     
%             IC0 = model.IC;
%     
%             model = update_grid(model, Parameters, Disp, allWfsn, truemodel);
%             model = inversion(model, Parameters, allWfsn, Disp, truemodel);
%     
%             dl = model.IC - IC0;
%     
%             disp([ ' --> Refined for an ' Parameters.IC ' of ' num2str(model.IC), ...
%                 '; improvement after z update of ' num2str(dl) ])
% 
%             if dl < 0
% 
%                 disp('     Fit can degrade when updating the z vector')
% 
%             end
%             
%             plot_model_new(5678, model, Parameters, allWfsn, Disp, truemodel);
%             plot_data(101, Parameters, model.synWfs, model.Disp);
%     
%             if dl < 1 && dl >= 0 %don't quit if its a z-update problem degrading the fit
% 
%                 keepgoing = false;
%     
%             end
%             
%             
%         end
%         
%     end
