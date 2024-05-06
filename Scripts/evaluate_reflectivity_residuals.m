function res = evaluate_reflectivity_residuals(model, allWfs, Disp, Parameters, errors, smth_flag)
    
    model_fine = fine_model(model, Parameters);
    %model_fine = fill_fields(model);

    lp = 0;

%     if any(model_fine.vs.model < Parameters.vs_lim(1))
% 
%         lp = -1e21;
% 
%     end
% 
%     if any(model_fine.vs.model > Parameters.vs_lim(2))
% 
%         lp = -1e22;
% 
%     end
% 
%     if any(model_fine.vpvs.model < Parameters.vpvs_lim(1))
% 
%         lp = -1e21;
% 
%     end
% 
%     if any(model_fine.vpvs.model > Parameters.vpvs_lim(2))
% 
%         lp = -1e22;
% 
%     end

    if model.sig_cr < Parameters.csig_limits(1) || model.sig_cr > Parameters.csig_limits(2)

        lp = -1e25;

    end

    if model.sig_cl < Parameters.csig_limits(1) || model.sig_cl > Parameters.csig_limits(2)

        lp = -1e25;

    end

    if any(model.sig_rf < Parameters.rf_limits(1)) || any(model.sig_rf > Parameters.rf_limits(2))

        lp = -1e26;

    end

    %if model.curvature < Parameters.curve_lim(1) || model.curvature > Parameters.curve_lim(2)

    %    lp = -1e30;

    %end

    %lp = lp - sum(0.5*(model.dvs.^2)/Parameters.vs_slope^2);

    if Parameters.debug || lp < -1e10

        llh  = 0;

    else

        [n,m] = size(allWfs);

        if Parameters.anirec && ~errors

            model.synWfs = do_rfs(model_fine, Parameters, allWfs);

        end

        ndata = 0;

        res = [];

        if Parameters.anirec
        
            for k = 1:length(Parameters.p_array)
    
                for j = 1:m
    
                    %includes filtering
    
                    %with weight
                    %phi_vec(k,j) = (allWfs(k,j).g/Parameters.sample_rate)*...
                    %    sum(((allWfs(k,j).rf - model.synWfs(k,j).rf)/(exp(model.sig_rf)/sqrt(model.synWfs(k,j).n))).^2);
                    %phi_vec(k,j) = model.synWfs(k,j).n*sum(((allWfs(k,j).rf - model.synWfs(k,j).rf)).^2);

                    ndata = ndata + 2*length(allWfs(k,j).rfr)*(allWfs(k,j).g/Parameters.sample_rate);

                    %n = n + model.synWfs(k,j).n;
                    res = [ res (allWfs(k,j).g/Parameters.sample_rate)*((allWfs(k,j).rfr - model.synWfs(k,j).rfr)')...
                        /(exp(model.sig_rf(j))/sqrt(model.synWfs(k,j).n)) ];
                    res = [ res (allWfs(k,j).g/Parameters.sample_rate)*((allWfs(k,j).rft - model.synWfs(k,j).rft)')...
                        /(exp(model.sig_rf(j))/sqrt(model.synWfs(k,j).n)) ];
     
                end
                
            end
    
        else
    
            phi_vecr = zeros(size(allWfs));
            phi_vect = zeros(size(allWfs));
            %res = 0;

        end

        %phi1 = sum(res)/n;%n weighted mean squared error

        %surface wave dispersion
        try

            if ~errors

                %model.cr_pred = dispR_surf96(Parameters.t_array, make_surf96model(model, Parameters));
                %model.cl_pred = zeros(size(model.cr_pred));
                %[model.cr_pred, model.cl_pred] = drive_aniprop(model_fine, Parameters);
                
                if ~any(model.vs.model < 0)
                
                    [model.cr_pred, model.cl_pred, ~, ~] = drive_mineos(model, Parameters, false);

                else 
                   
                    model.cr_pred = 1e9*ones(size(Parameters.t_array));
                    model.cl_pred = 1e9*ones(size(Parameters.t_array));

                end

                if isempty(model.cr_pred)

                    model.cr_pred = 1e9*ones(size(Parameters.t_array));
                    model.cl_pred = 1e9*ones(size(Parameters.t_array));

                end



            end

            if ndata > 0 || ~Parameters.anirec

                n = length(Disp.c_r) + length(Disp.c_l);

            else

                n     = 1;
                ndata = 4;
                
            end

            %greatly upweight the surface wave data
            res   = [ res 0.25*(ndata/n)*((model.cr_pred - Disp.c_r)/exp(model.sig_cr))' ];
            res   = [ res 0.25*(ndata/n)*((model.cl_pred - Disp.c_l)/exp(model.sig_cl))' ];
            %phi2   = sum((model.cr_pred - Disp.c_r).^2)/length(Disp.c_r);
            %phi2   = sum(((model.cr_pred - Disp.c_r)/exp(model.sig_cr)).^2);
            %llh    = -0.5*(phi1/exp(model.sig_rf)^2 + phi2/exp(model.sig_cr)^2);
            %llh    = -0.5*(sum(phi_vec(:)) + phi2);

        catch

            res   = [ res 1e9*ones(size(Disp.c_r))' ];
            res   = [ res 1e9*ones(size(Disp.c_l))' ];
            %llh  = -1e20;

        end

        z       = (0:0.1:Parameters.plot_z)';
        vs_sm   = interp1(model.vs.z,   model.vs.node, z, 'linear', 'extrap');
        vpvs_sm = interp1(model.vpvs.z, model.vpvs.node, z, 'linear', 'extrap');

        vs_sm(z<model.vs.z(1))     = model.vs.node(1);
        vpvs_sm(z<model.vpvs.z(1)) = model.vpvs.node(1);

        res = [ res ((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2)) ];

        %smoothness?
        if smth_flag
                
            res = [ res gradient((vs_sm' - Parameters.limits.vs(1))/diff(Parameters.limits.vs), 0.1)/exp(model.curvature) ];        
            res = [ res gradient((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), 0.1)/exp(model.curvature) ];

        end

    end
    
end
