function model = reduce_grid(model, Parameters, Disp, allWfs, truemodel)

    keepgoing = true;

    while keepgoing

        keepgoing = false;

        for fidx = 1:length(Parameters.fields_vec)
    
            field = Parameters.fields_vec{fidx};
    
            %remove nodes from the bottom of the mesh first
            k = model.(field).n-1;
    
            disp([ 'Reducing the ' field ' grid....' ])
            
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
    
                %mn = update_errors(mn, allWfs, Disp, Parameters);

                %mn = update_z(mn, Parameters);
                mn = evaluate_reflectivity(mn, allWfs, Disp, Parameters, false);
    
                if (model.IC - mn.IC) < 0
    
                    disp([ '     Node removed at ' num2str(model.(field).z(k)*model.max_z) ' km; ' Parameters.IC ' now ' num2str(mn.IC) ])
                    %mn    = update_z(mn, Parameters);
                    %mn    = evaluate_reflectivity(mn, allWfs, Disp, Parameters, false);
                    model = mn;
                    plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
                    %plot_data(101, Parameters, model.synWfs, model.Disp);
    
                    keepgoing = true;

                end
    
                k = k-1;
    
            end
    
        end

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

