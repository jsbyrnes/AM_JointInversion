function model = update_grid(model, Parameters, Disp, allWfs, truemodel)

    for fidx = 1:length(Parameters.fields_vec)

        field = Parameters.fields_vec{fidx};

        %remove nodes from the mesh first
        k = model.(field).n-1;

        disp([ 'Updating the ' field ' grid....' ])
        
        while k > 1 && model.(field).n > 3

            mn = model;

            mn.(field).node(k) = [];
            mn.(field).z(k)    = [];
            mn.(field).n       = mn.(field).n - 1;

            dz_check = max([ 0; diff(mn.(field).z)*model.max_z ]);

            if dz_check > Parameters.dz_lim(2)

                k = k-1;
                continue

            end

            mn = update_z(mn, Parameters);
            mn = evaluate_reflectivity(mn, allWfs, Disp, Parameters, false);

            if (model.IC - mn.IC) < 0

                disp([ '     Node removed at ' num2str(model.(field).z(k)*model.max_z) ' km; ' Parameters.IC ' now ' num2str(mn.IC) ])
                mn    = update_z(mn, Parameters);
                mn    = evaluate_reflectivity(mn, allWfs, Disp, Parameters, false);
                model = mn;
                plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
                plot_data(101, Parameters, model.synWfs, model.Disp);

            end

            k = k-1;

        end

        %add them back in with small 1-d optimization problems
        %first make a list of all the new possible depths
        zlist = [];

        for k = 1:(model.(field).n-1)

            if k == 1 && model.(field).z(1)*model.max_z > 0 %need to use if the top node isn't the surface

                dz = model.(field).z(1)*model.max_z;

            else

                dz = (model.(field).z(k+1) - model.(field).z(k))*model.max_z;

            end

            zlist = [ zlist; (model.(field).z(k)*model.max_z + ((1:Parameters.division-1)')*dz/Parameters.division) ];

        end

        zlist(zlist<=0)                 = [];
        zlist(zlist>=Parameters.plot_z) = [];
        zlist                           = unique(sort(zlist));
        model                           = fine_model(model, Parameters);

        for k = 1:length(zlist)

            mn = model;

            mn.(field).z    = [ mn.(field).z; zlist(k)/model.max_z ];
            mn.(field).node = [ mn.(field).node; interp1(model.z, model.(field).model, zlist(k), Parameters.interp, 'extrap') ];

            [mn.(field).z, ix] = sort(mn.(field).z);
            mn.(field).node = mn.(field).node(ix);

            mn = update_z(mn, Parameters);

            %and find it
            [~, ix] = min(abs(mn.(field).z*model.max_z - zlist(k)));
    
            %is it too close to stuff?
            dz_check = min(abs(model.(field).z*model.max_z - zlist(k)));
            if dz_check < Parameters.dz_lim(1)

                continue

            end

            func = @(x) wrapper(x, mn, Disp, allWfs, Parameters, field, ix);
            
            dx = 1;
            x0 = mn.(field).node(ix);
            f0 = func(x0);

            %iterate newton's method (maximization)
            count = 0;
            while abs(dx) > 1e-2 && count < 4 && f0 < model.IC

                fp = func(x0 + Parameters.h);
                fn = func(x0 - Parameters.h);

                d1 = (0.5*fp - 0.5*fn)/Parameters.h;
                d2 = (fp + fn - 2*f0)/Parameters.h^2;

                dx = d1/d2;

                if abs(dx) > 0.5

                    dx = 0.5*sign(d1);

                end

                x0 = x0 - dx;%maximization

                %enforce bounds
                if strcmp(field, 'vs')

                   x0 = min([ max([ x0 Parameters.limits.vs(1) ]) Parameters.limits.vs(2) ]); 

                elseif strcmp(field, 'vpvs')

                   x0 = min([ max([ x0 1 ]) 5 ]);%models will start to fail if it gets unrealistic

                end                
                
                f0    = func(x0);
                count = count + 1;

            end
            
            if f0 - model.IC > 0

                mn.(field).node(ix) = x0;
                mn                  = update_z(mn, Parameters);
                mn                  = evaluate_reflectivity(mn, allWfs, Disp, Parameters, false);
                model               = mn;
                model               = fine_model(model, Parameters);

                disp([ '     Node added at ' num2str(mn.(field).z(ix)*model.max_z) ' km; ' Parameters.IC ' now ' num2str(mn.IC) ])
                plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
                plot_data(101, Parameters, model.synWfs, model.Disp);
                
            end

        end

        %model.(field).z = model.(field).z/model.max_z;

    end

end

%     %update the errors
%     func = @(x) wrapper_sig(x, model, Disp, allWfs, Parameters);
%     opt  = optimset('MaxIter', 1e100, 'MaxFunEvals', 1e100);
%     vec  = fminsearch(func, model.sig_rf, opt);%NM handles bounds well
% 
%     model.sig_rf = vec;
%     model        = evaluate_reflectivity(model, allWfs, Disp, Parameters, true);

%             x0   = mn.(field).node(ix);
% 
%             xsearch = linspace(-0.5, 0.5, Parameters.numworkers);
% 
%             parfor j = 1:length(xsearch)
% 
%                 fsearch(j) = func(x0 + xsearch(j));
% 
%             end
% 
%             xn = -0.5:0.001:0.5;
%             fn = interp1(xsearch, fsearch, xn, 'spline');
%             [~, fix] = min(fn);
%             x0 = x0 + xn(fix);
%             f0 = func(x0);
% 

            %[ mn.(field).node(ix), fval ] = fminbnd(func, mn.(field).node(ix)-0.25, mn.(field).node(ix)+0.25, opt_single);

    %update the errors
%     func = @(x) wrapper_sig(x, model, Disp, allWfs, Parameters);
%     vec  = fminsearch(func, [ model.sig_cr model.sig_cl model.sig_rf ]);%NM handles bounds well
% 
%     model.sig_cr = vec(1);
%     model.sig_cl = vec(2);
%     model.sig_rf = vec(3:(3+length(Parameters.g_array)-1));
%     model        = evaluate_reflectivity(model, allWfs, Disp, Parameters, true);

%     y = -6:0.01:0; 
%     for k = 1:length(y)
%         
%         llh(k) = -y(k)*length(Disp.c_r) - 0.5*sum(((model.cr_pred - Disp.c_r).^2)/exp(y(k))^2);
% 
%     end
% 
%     [~, ix] = max(llh);
%     model.sig_cr = y(ix);
% 
%     for k = 1:length(y)
%         
%         llh(k) = -y(k)*length(Disp.c_l) - 0.5*sum(((model.cl_pred - Disp.c_l).^2)/exp(y(k))^2);
% 
%     end
% 
%     [~, ix] = max(llh);
%     model.sig_cl = y(ix);
% 
%     %set it to the std
%     for g = 1:length(Parameters.g_array)
% 
%         res = [];
% 
%         for k = 1:length(Parameters.p_array)
%             
%             res = [ res (allWfs(k,g).rfr - model.synWfs(k,g).rfr)' ];
%             res = [ res (allWfs(k,g).rft - model.synWfs(k,g).rft)' ];
% 
%         end
% 
%         model.sig_rf(g) = log(std(res));
% 
%     end
% 
%     model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false); %update the data

