function model = evaluate_reflectivity(model, allWfs, Disp, Parameters, errors)

    model = fine_model(model, Parameters);

    for k = 1:length(Parameters.fields_vec)

        model.(Parameters.fields_vec{k}).n   = length(model.(Parameters.fields_vec{k}).node);

        %dimensionalize the depths
        model.(Parameters.fields_vec{k}).z = model.(Parameters.fields_vec{k}).z*model.max_z;

    end

    lp = 0;

%     if any(model_fine.vs.model < Parameters.limits.vs(1))
% 
%         model_fine.vs.model(model_fine.vs.model < Parameters.limits.vs(1)) = Parameters.limits.vs(1);
% 
%     end
% 
%     if any(model_fine.vs.model > Parameters.limits.vs(2))
% 
%         model_fine.vs.model(model_fine.vs.model > Parameters.limits.vs(1)) = Parameters.limits.vs(2);
% 
%     end

    if isfield(model, 'parameterized')

        lp = lp - 0.5*sum(model.parameterized.vec.^2);%already normalized
        lp = lp - 0.5*(( (model.parameterized.mantle(1) - model.parameterized.crust(end)) - 0.5 )^2)/0.2^2;

    end

    if Parameters.debug || lp < -1e10

        llh  = 0;

    else

        if Parameters.anirec && ~errors

            model.synWfs = do_rfs(model, Parameters, allWfs);

        end

        ndata    = 0;%zeros(size(Parameters.g_array));
        rffactor = 0;%zeros(size(Parameters.g_array));

        res      = [];

        if isfield(model, 'synWfs')
        
            for k = 1:length(Parameters.p_array)
    
                for j = 1:length(Parameters.g_array)

                    phi_vec(k,j) = (allWfs(k,j).g/Parameters.sample_rate)*sum(((allWfs(k,j).rfr - model.synWfs(k,j).rfr)./...
                        allWfs(k,j).rfr_std).^2);

                    %rffactor = rffactor - length(allWfs(k,j).rfr)*(allWfs(k,j).g/Parameters.sample_rate)...
                    %    *log(exp(model.sig_rf(j))/sqrt(model.synWfs(k,j).n));%/sqrt(model.synWfs(k,j).n)

                    rffactor = rffactor - (allWfs(k,j).g/Parameters.sample_rate)...
                        *sum(log(allWfs(k,j).rfr_std));

                    res = [ res sqrt(allWfs(k,j).g/Parameters.sample_rate)*(((allWfs(k,j).rfr - model.synWfs(k,j).rfr))...
                        ./allWfs(k,j).rfr_std)' ];
    
                end
                
            end
    
        else
    
            phi_vec = zeros(size(allWfs));
            %res = 0;

        end

        %phi1 = sum(res)/n;%n weighted mean squared error

        %surface wave dispersion
        if ~errors && Parameters.surf_96

            try

                [ cr_pred, t_found ] = dispR_surf96(Parameters.t_array, make_surf96model(model, Parameters));

                %in case there is a mode jump. Will only occur if model is
                %getting discarded anyway
                %99% of the time this does nothing
                model.cr_pred = interp1(t_found, cr_pred, Parameters.t_array, 'nearest', 'extrap');

                if isrow(model.cr_pred)

                    model.cr_pred = model.cr_pred';

                end

            catch

                model.cr_pred = model.cr_pred;%dummy statement

            end

            %model.cl_pred = zeros(size(model.cr_pred));
            %[model.cr_pred, model.cl_pred] = drive_aniprop(model, Parameters);
            
            %if ~any(model.vs.model < 0)
            %
            %    [model.cr_pred, model.cl_pred, ~, ~] = drive_mineos(model, Parameters, false);
            %
            %else 
            %  
            %    model.cr_pred = 1e9*ones(size(Parameters.t_array));
            %    model.cl_pred = 1e9*ones(size(Parameters.t_array));
            %
            %end
           
            if isempty(model.cr_pred)

                model.cr_pred = 1e9*ones(length(Parameters.t_array), 1);
                model.cl_pred = 1e9*ones(length(Parameters.t_array), 1);

            end

        end

        phi2   = sum(((model.cr_pred - Disp.c_r)./Disp.c_rstd).^2);% + sum(((model.cl_pred - Disp.c_l)/exp(model.sig_cl)).^2);
        llh    = -0.5*(sum(phi_vec(:)) + phi2);

        if any(ndata > 0) && Parameters.anirec

            n = length(Disp.c_r);% + length(Disp.c_l);

        else

            n     = 1;
            ndata = 4;
            
        end

        %res   = [ res ((Disp.c_r - model.cr_pred)./Disp.c_rstd)' ];
        res   = [ res ((Disp.c_r - model.cr_pred)./Disp.c_rstd)' ];%
           
        %llh = llh - sum(model.sig_rf(1:length(ndata)).*ndata) - model.sig_cr*length(Disp.c_r);% - model.sig_cl*length(Disp.c_l);% - model.sig_cl*length(Disp.c_l); 
        %llh = llh - model.sig_rf*sum(ndata);% - model.sig_cl*length(Disp.c_l);% - model.sig_cl*length(Disp.c_l); 
        llh  = llh - sum(log(Disp.c_rstd)) + rffactor;
        %llh = llh - model.sig_rf - model.sig_cr;

        z       = (0:Parameters.dz_lim(1):Parameters.plot_z)';
        %vs_sm   = interp1(model.vs.z,   model.vs.node, z, 'linear', 'extrap');
        %vpvs_sm = interp1(model.vpvs.z, model.vpvs.node, z, 'linear', 'extrap');
        vs_sm   = interp1(model.vs.z, model.vs.node, z, Parameters.interp, model.vs.model(end));

        %smoothness        
        if ~isempty(Parameters.Firstgrad_damp)
        
            %lp = lp - 0.5*sum((gradient((vs_sm - Parameters.limits.vs(1))/diff(Parameters.limits.vs), Parameters.dz_lim(1))/exp(Parameters.Firstgrad_damp)).^2);        
            %res = [ res -1*gradient((vs_sm' - Parameters.limits.vs(1))/diff(Parameters.limits.vs), Parameters.dz_lim(1))/exp(Parameters.Firstgrad_damp) ];        

                vs       = exp(model.vs.node);
                z        = model.vs.z;
                grad1     = diff(vs)./diff(z);

                %weaker at the top, but punished if positive at the top
                w = exp(Parameters.Firstgrad_damp)*ones(size(grad1)).*(1 - 0.99*exp(-model.vs.z(1:length(grad1))/(Parameters.sed_model_z)));

                %ix          = model.vs.z(2:end) <= Parameters.sed_model_z;
                %grad(ix)    = grad(ix)/10;
                %continuous with the parameterized top
                %ix1      = find(model.z<Parameters.sed_model_z, 1, 'last');
                %ix2      = find(model.z>Parameters.sed_model_z, 1, 'first');
                %topgrad  = (model.vs.model(ix2) - model.vs.model(ix1))/(model.z(ix2) - model.z(ix1));
                %grad     = diff(diff([ model.vs.model(ix1-1); model.vs.model(ix1); vs])./diff([ model.z(ix1-1); model.z(ix1); z ]));
                %grad(1) = grad(1)*100;

                %ix = model.vs.z(2:end) <= Parameters.sed_model_z;
                %grad(ix) = 0;
                    
                %grad_weight = ones(length(grad),1);
                %grad_weight(model.vs.z(2:(length(grad)+1)) < Parameters.plot_z)  = exp(Parameters.Secondgrad_damp);
                %grad_weight(model.vs.z(2:(length(grad)+1)) >= Parameters.plot_z) = 5*exp(Parameters.Secondgrad_damp);

                if strcmp(Parameters.smoothing_style, 'cauchy')

                    lp      = lp - 0.5*sum(log(1 + (grad1.*w).^2));
                    res     = [ res -1*sqrt(abs(log(1 + (grad1.*w).^2)))' ];

                elseif strcmp(Parameters.smoothing_style, 'laplacian')

                    lp      = lp - 0.5*sum(abs(grad1.*w));
                    res     = [ res -1*sqrt(abs((grad1.*w)))' ];

                elseif strcmp(Parameters.smoothing_style, 'gaussian')

                    lp      = lp - 0.5*sum((grad1.*w).^2);
                    res     = [ res -1*((grad1.*w))' ];
    
                end

        end

        if ~isempty(Parameters.Secondgrad_damp)

            if isfield(model, 'parameterized')
        
                %sediment gradients
                vs      = model.vs.node(model.vs.z <= exp(model.parameterized.sed_z));
                z       = model.vs.z(model.vs.z <= exp(model.parameterized.sed_z));
                grad    = diff(diff(vs)./diff(z));
                lp      = lp - 0.5*sum((grad*exp(Parameters.Secondgrad_damp)).^2);
                res     = [ res -1*(grad*exp(Parameters.Secondgrad_damp)) ];        

                %crust gradients
                vs      = model.vs.node(model.vs.z >= exp(model.parameterized.sed_z) & model.vs.z <= exp(model.parameterized.moho_z));
                z       = model.vs.z(model.vs.z >= exp(model.parameterized.sed_z) & model.vs.z <= exp(model.parameterized.moho_z));
                grad    = diff(diff(vs)./diff(z));
                lp      = lp - 0.5*sum((grad*exp(Parameters.Secondgrad_damp)).^2);
                res     = [ res -1*(grad*exp(Parameters.Secondgrad_damp)) ];        

                %mantle gradients
                vs      = model.vs.node(model.vs.z >= exp(model.parameterized.moho_z) + exp(model.parameterized.moho_dz));
                z       = model.vs.z(model.vs.z >= exp(model.parameterized.moho_z) + exp(model.parameterized.moho_dz));
                grad    = diff(diff(vs)./diff(z));
                lp      = lp - 0.5*sum((grad*exp(Parameters.Secondgrad_damp)).^2);
                res     = [ res -1*(grad*exp(Parameters.Secondgrad_damp)) ];        

            else

                %lp  = lp - 0.5*sum(log(1 + (gradient(gradient((vs_sm - Parameters.limits.vs(1))...
                %   /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp)).^2));
                
                %res = [ res -1*sqrt(abs(log(1 + (gradient(gradient((vs_sm - Parameters.limits.vs(1))...
                %   /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp)).^2)))' ];

                %lp  = lp - 0.5*sum(abs(gradient(gradient((vs_sm - Parameters.limits.vs(1))...
                %    /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp)));
                
                %res = [ res -1*sqrt(abs(gradient(gradient((vs_sm' - Parameters.limits.vs(1))...
                %    /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp))) ];

                %lp  = lp - 0.5*sum((gradient(gradient((vs_sm - Parameters.limits.vs(1))...
                %   /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp)).^2);
                
                %res = [ res -1*gradient(gradient((vs_sm' - Parameters.limits.vs(1))...
                %   /diff(Parameters.limits.vs), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp) ];

                vs       = exp(model.vs.node);
                z        = model.vs.z;
                grad     = diff(vs)./diff(z);
                grad2    = diff(grad);

                %weaker at the top, but punished if positive at the top
                w = exp(Parameters.Secondgrad_damp)*ones(size(grad2)).*(1 - 0.99*exp(-model.vs.z(1:length(grad2))/(Parameters.sed_model_z)));

                %ix          = model.vs.z(2:end) <= Parameters.sed_model_z;
                %grad(ix)    = grad(ix)/10;
                %continuous with the parameterized top
                %ix1      = find(model.z<Parameters.sed_model_z, 1, 'last');
                %ix2      = find(model.z>Parameters.sed_model_z, 1, 'first');
                %topgrad  = (model.vs.model(ix2) - model.vs.model(ix1))/(model.z(ix2) - model.z(ix1));
                %grad     = diff(diff([ model.vs.model(ix1-1); model.vs.model(ix1); vs])./diff([ model.z(ix1-1); model.z(ix1); z ]));
                %grad(1) = grad(1)*100;

                %ix = model.vs.z(2:end) <= Parameters.sed_model_z;
                %grad(ix) = 0;
                    
                %grad_weight = ones(length(grad),1);
                %grad_weight(model.vs.z(2:(length(grad)+1)) < Parameters.plot_z)  = exp(Parameters.Secondgrad_damp);
                %grad_weight(model.vs.z(2:(length(grad)+1)) >= Parameters.plot_z) = 5*exp(Parameters.Secondgrad_damp);

                if strcmp(Parameters.smoothing_style, 'cauchy')

                    lp      = lp - 0.5*sum(log(1 + (grad2.*w).^2));
                    res     = [ res -1*sqrt(abs(log(1 + (grad2.*w).^2)))' ];

                elseif strcmp(Parameters.smoothing_style, 'laplacian')

                    lp      = lp - 0.5*sum(abs(grad2*exp(Parameters.Secondgrad_damp)));
                    res     = [ res -1*sqrt(abs((grad2*exp(Parameters.Secondgrad_damp))))' ];

                elseif strcmp(Parameters.smoothing_style, 'gaussian')

                    lp      = lp - 0.5*sum((grad2*exp(Parameters.Secondgrad_damp)).^2);
                    res     = [ res -1*((grad2*exp(Parameters.Secondgrad_damp)))' ];
    
                end

            end

        end

        if any(contains(Parameters.fields_vec, 'vpvs'))

            %vpvs_sm = interp1(model.z, model.vpvs.model, z, Parameters.interp, model.vpvs.model(1));
            lp  = lp - 0.5*sum(((model.vpvs.node - 1.76)/0.1).^2);
            res = [ res -1*((model.vpvs.node - 1.76)/0.1)' ];

            vpvs    = model.vpvs.node;
            z       = model.vpvs.z;
            grad    = diff(diff(vpvs)./diff(z));

            %lp      = lp - 0.5*sum(log(1 + (grad*exp(Parameters.Secondgrad_damp)).^2));
            %res     = [ res -1*sqrt(abs(log(1 + (grad*exp(Parameters.Secondgrad_damp)).^2)))' ];

%             %smoothness        
%             if ~isempty(Parameters.Firstgrad_damp)
%             
%                 lp  = lp - 0.5*sum((gradient((vpvs_sm - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), Parameters.dz_lim(1))/exp(Parameters.Firstgrad_damp)).^2);
%                 res = [ res -1*gradient((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), Parameters.dz_lim(1))/exp(Parameters.Firstgrad_damp) ];
%     
%             end
%     
%             if ~isempty(Parameters.Secondgrad_damp)
%     
%                 lp  = lp - 0.5*sum((gradient(gradient((vpvs_sm - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp)).^2);
%                 res = [ res -1*gradient(gradient((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), Parameters.dz_lim(1)), Parameters.dz_lim(1))*exp(Parameters.Secondgrad_damp) ];
% 
%             end
    
        end

        %
        %and surface and depth contraints
        %lp  = lp - 0.5*((model.vs.model(1) - model.vs_surface(1))^2)/model.vs_surface(2)^2;
        %res = [ res (model.vs.model(1) - model.vs_surface(1))/model.vs_surface(2) ];

        %%%%%%%%%%%%%%
        %monotonicity constraint on the top
        grad = diff(exp(model.vs.node));
        dz   = diff(model.vs.z);
        sed_grads   = (grad./dz)';
        sed_grads(sed_grads < -20) = -20;%very rare but avoids and inifinty
        sed_lp      = exp(-sed_grads*Parameters.monotonicity_w).*exp(-model.vs.z(1:length(sed_grads))/(Parameters.sed_model_z))';

        lp  = lp - 0.5*sum(sed_lp.^2);
        res = [ res -sed_lp ];
        %%%%%%%%%%%%%%

        if model.vpvs_on

            lp  = lp - 0.5*((model.vpvs_z - log(0.25))/0.25)^2;
            res = [ res ((model.vpvs_z - log(0.25))/0.25) ];

            lp  = lp - 0.5*((model.vpvs_block - log(1.76))/0.25)^2;
            res = [ res ((model.vpvs_block - log(1.76))/0.25) ];

        end

        %priors/res on sed model
%         lp = lp - 0.5*((model.sed_model(1) - log(0.25))/0.5)^2;
%         lp = lp - 0.5*((model.sed_model(2) - log(2))/0.75)^2;
%         lp = lp - 0.5*((model.sed_model(3) - log(2))/2.5)^2;
%         lp = lp - 0.5*((model.sed_model(4) - log(1))/1)^2;
%         lp = lp - 0.5*(model.sed_model(5) - log(0.2/Parameters.sed_model_z))^2/(0.5^2);
% 
%         res = [ res (model.sed_model(1) - log(0.25))/0.5 (model.sed_model(2) - log(2))/0.75 (model.sed_model(3) - log(2))/2.5 ...
%             (model.sed_model(4))/1 (model.sed_model(5) - log(0.2/Parameters.sed_model_z))/0.5 ];

        lp  = lp - 0.5*((model.max_z - Parameters.plot_z)/Parameters.zstd)^2;
        res = [ res ((model.max_z - Parameters.plot_z)/Parameters.zstd) ];

%         lp  = lp - 0.5*((model.vs.model(1) - model.vs0)/0.01)^2;
%         res = [ res ((model.vs.model(1) - model.vs0)/0.01) ];

    end

    model.res    = res;
    model.L2norm = sum(res.^2);
    model.lpst   = llh + lp;
    model.lp     = lp;
    model.llh    = llh;
    model.n      = sum(ndata) + length(Disp.c_r);

    %redundent for plotting
    model.Disp     = Disp;
    model.Disp.c_r = model.cr_pred;

    n = 0;
    
    for k = 1:length(Parameters.fields_vec)

        n = n + model.(Parameters.fields_vec{k}).n;

    end

    if model.vpvs_on

        n = n + 2;

    end

    if strcmp(Parameters.IC, 'AIC')

        model.IC = model.lpst - n*2;%log(sum(ndata) + length(Disp.c_r));

    elseif strcmp(Parameters.IC, 'BIC')

        model.IC = model.lpst - n*log(sum(ndata) + length(Disp.c_r));

    elseif strcmp(Parameters.IC, 'lnIC')

        model.IC = model.lpst - log(n);

    end

    for k = 1:length(Parameters.fields_vec)

        %re-nondimensionalize the depths
        model.(Parameters.fields_vec{k}).z = model.(Parameters.fields_vec{k}).z/model.max_z;

    end

end

%         vs   = (model.vs.node - Parameters.limits.vs(1))/diff(Parameters.limits.vs);
%         vpvs = ((model.vpvs.node - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2));
% 
%         lp = lp - 0.5*sum(vpvs.^2);
%         res = [ res vpvs' ];
% 
%         ddvsddz   = (diff([vs(1); vs])./diff([ 0; model.vs.z]));
%         %ddvpvsddz = diff(diff([vpvs(1); vpvs])./diff([ 0; model.vpvs.z]));
%         dvpvsdz = diff([vpvs(1); vpvs])./diff([ 0; model.vpvs.z]);
% 
%         lp  = lp - 0.5*sum(ddvsddz.^2)/exp(model.curvature)^2;
%         %lp  = lp - 0.5*sum(ddvpvsddz.^2)/exp(model.curvature)^2;
%         lp  = lp - 0.5*sum(dvpvsdz.^2)/0.01^2;
%         res = [ res (ddvsddz/exp(model.curvature))' ];
%         %res = [ res (ddvpvsddz/exp(model.curvature))' ];
%         res = [ res (dvpvsdz/0.01)' ];
%         res = [ res ((model.vs.model(1) - model.vs_surface(1)))/model.vs_surface(2) ];

%         z       = (0:0.1:Parameters.plot_z)';
%         %vs_sm   = interp1(model.vs.z,   model.vs.node, z, 'linear', 'extrap');
%         %vpvs_sm = interp1(model.vpvs.z, model.vpvs.node, z, 'linear', 'extrap');
%         vs_sm   = interp1(model.z, model.vs.model, z, Parameters.interp, model.vs.model(1));
%         vpvs_sm = interp1(model.z, model.vpvs.model, z, Parameters.interp, model.vpvs.model(1));
% 
%         %damp
%         lp = lp - 0.5*sum(((vpvs_sm - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2)).^2);
% 
%         %smoothness
%         %lp = lp - 0.5*sum((gradient(gradient((vs_sm - Parameters.limits.vs(1))/diff(Parameters.limits.vs), 0.1), 0.1)/exp(model.curvature)).^2);        
%         %lp = lp - 0.5*sum((gradient(gradient((vpvs_sm - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), 0.1), 0.1)/exp(model.curvature)).^2);
%         lp = lp - 0.5*sum((gradient((vs_sm - Parameters.limits.vs(1))/diff(Parameters.limits.vs), 0.1)/exp(model.curvature)).^2);        
%         lp = lp - 0.5*sum((gradient((vpvs_sm - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), 0.1)/exp(model.curvature)).^2);
% 
%         %and surface constraint
%         lp = lp - 0.5*((model.vs.model(1) - model.vs_surface(1))^2)/model.vs_surface(2)^2;
% 
%         %repeat as res
%         %vs_sm(z<model.vs.z(1))     = model.vs.node(1);
%         %vpvs_sm(z<model.vpvs.z(1)) = model.vpvs.node(1);
% 
%         res = [ res ((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2)) ];
% 
%         %res = [ res gradient(gradient((vs_sm' - Parameters.limits.vs(1))/diff(Parameters.limits.vs), 0.1), 0.1)/exp(model.curvature) ];        
%         %res = [ res gradient(gradient((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), 0.1), 0.1)/exp(model.curvature) ];
%         res = [ res gradient((vs_sm' - Parameters.limits.vs(1))/diff(Parameters.limits.vs), 0.1)/exp(model.curvature) ];        
%         res = [ res gradient((vpvs_sm' - Parameters.limits.vpvs(1))/Parameters.limits.vpvs(2), 0.1)/exp(model.curvature) ];


%%%%%%%%
%old search for alignment
%                     lags = 1:length(dN);
%                     for j = 1:length(lags)
% 
%                         %comutationally intensive
%                         phi(j) = md(circshift(No, lags(j)), circshift(Eo, lags(j)));
%     
%                     end
%     
%                     [R(2), idx] = min(phi);
%                     %find neighbors, assuming periodic
%                     if idx==1
%     
%                         R(1)=phi(end);
%                         R(3)=phi(2);
%     
%                     elseif idx==numel(phi)
%     
%                         R(1)=phi(end-1);
%                         R(3)=phi(1);
%     
%                     else
%     
%                         R(1)=phi(idx-1);
%                         R(3)=phi(idx+1);
%     
%                     end
%     
%                     c     = (R(3)-R(1))/(2*(2*R(2)-R(1)-R(3)));
%                     %lag   = mod(idx-1+floor(length(lags)/2),length(lags))-floor(length(lags)/2);%integer part
%                     delay = idx+c;%delay estimate
%                     N     = delay_continuous(No, Parameters.sample_rate, delay/Parameters.sample_rate );
%                     E     = delay_continuous(Eo, Parameters.sample_rate, delay/Parameters.sample_rate );


%double grid search for amplitude
%                     phi = []
%                     amps = (-6:0.05:3)
%                     for j = 1:length(amps)
%     
%                         phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%     
%                     end
%     
%                     [mphi, idx] = min(phi);
% 
%                     if idx > 1
% 
%                         %refined search - needs to be extremely accurate
%                         phi = [];
%                         amps = (-0.5:0.05:0.5) + amps(idx);
%                         for j = 1:length(amps)
%         
%                             %computationally intensive
%                             phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%         
%                         end
%                         [mphi, idx] = min(phi);        
%                         c     = (phi(idx+1) - phi(idx-1))/(2*(2*mphi - phi(idx-1) - phi(idx+1)));
%                         %ampi  = mod(idx - 1 + floor(length(amps)/2), length(amps)) - floor(length(amps)/2);%integer part
%                         amp   = amps(idx) + c*0.05;%amp estimate
% 
%                     else
% 
%                         amp = -6;%basically zero, its a huge mismatch. 
% 
%                     end
% 
%                     N = exp(amp)*N;
%                     E = exp(amp)*E;

%         if sum(model.dt(k, :)) > model.dt(k, 1)%check for a second layer
% 
%             %account for rotation
% %             dphi     = abs( (model.fast_dir(k, 1) + model.fast_dir_rotation(k,1)) ...
% %                 - (model.fast_dir(k, 2) - model.fast_dir_rotation(k,2)));
%             dphi     = abs( model.fast_dir(k, 1) - model.fast_dir(k, 2));
% 
%             if dphi > pi/2
% 
%                 dphi = dphi - pi/2;
% 
%             end
% 
%             dphi = dphi/(pi/2);
% 
%             model.lp = model.lp + (Parameters.beta - 1)*(log(dphi) + log(1 - dphi));
% 
%         end

%     if nsta > 1
% 
%         for k = 1:length(model.sources)
%     
%             %model.lp = model.lp - 0.5*sum(model.sources{k}.^2);
%     
%     %         model.lp = model.lp - 0.5*sum(model.sA{k}.^2);
%     %         model.lp = model.lp - 0.5*sum(model.sB{k}.^2);
%     % 
%     %         model.lp = model.lp - sum(exp(-model.st{k}));
%     %         model.lp = model.lp - sum(exp((model.st{k} - model.t(end))));
%     %         model.lp = model.lp - sum(exp(-(model.sf{k} - minf)));
%     %         model.lp = model.lp - sum(exp((model.sf{k} - maxf)));
%     
%             ns = ns + length(model.sources{k});
%     
%         end
% 
%     end

%     if Parameters.use_covarience
% 
%         model.lp = model.lp - (((model.r - Parameters.r_range(1))^2)/Parameters.r_range(2)^2 + ...
%             sum(((model.f - Parameters.f_range(1)).^2)./Parameters.f_range(2)^2))/2;
%         %model.lp = model.lp - 1e9*sum(model.a<0);
%         %model.lp = model.lp - 1e9*sum(model.a>1);
% 
%         if exp(model.f) > Parameters.low_pass
%             
%             model.lp = -1e9;
%     
%         end
% 
%     end

%     if model.sig_cr < Parameters.csig_limits(1) || model.sig_cr > Parameters.csig_limits(2)
% 
%         lp = -1e25;
% 
%     end
% 
%     %if model.sig_cl < Parameters.csig_limits(1) || model.sig_cl > Parameters.csig_limits(2)
%     %
%     %    lp = -1e25;
%     %
%     %end
% 
%     if any(model.sig_rf < Parameters.rf_limits(1)) || any(model.sig_rf > Parameters.rf_limits(2))
% 
%         lp = -1e26;
% 
%     end
% 
%     if model.curvature < Parameters.curve_lim(1) || model.curvature > Parameters.curve_lim(2)
% 
%         lp = -1e30;
% 
%     end
% 
