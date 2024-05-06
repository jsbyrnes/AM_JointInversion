function model = expand_grid(model, Parameters, Disp, allWfs, probability)

    if nargin == 4

        probability = 1;%never skip

    end

    model = fine_model(model, Parameters);

    %will break if you add in other fields
    for fidx = 1:length(Parameters.fields_vec)

        field = Parameters.fields_vec{fidx};

        zlist = [];

        %model.(field).z    = [ Parameters.sed_model_z/Parameters.plot_z; model.vs.z ];    
        %model.(field).node = [ interp1(model.z/Parameters.plot_z, model.(field).model, model.(field).z(1)); model.(field).node ];

        %seperately do the top
        dz    = model.(field).z(1)*model.max_z;

        if dz/Parameters.division > Parameters.dz_lim(1)

            zlist = [ zlist; ((1:Parameters.division-1)')*dz/Parameters.division ];

        end

        %treat these differently
        upper_z  = model.(field).z(1)*model.max_z;
        upper_vs = exp(model.(field).node(1));

        for k = 2:(model.(field).n)

%             if k == 1 && model.(field).z(1)*Parameters.plot_z > 0 %need to use if the top node isn't the surface
% 
%                 dz = model.(field).z(1)*Parameters.plot_z;
% 
%             else
% 
            dz = (model.(field).z(k) - model.(field).z(k-1))*model.max_z;
% 
%             end

            if dz/Parameters.division > Parameters.dz_lim(1)

                zlist = [ zlist; (model.(field).z(k-1)*model.max_z + ((1:Parameters.division-1)')*dz/Parameters.division) ];

            end

        end

        %model.(field).z    = model.(field).z(2:end);
        %model.(field).node = model.(field).node(2:end);

        zlist(zlist<=0)                 = [];
        zlist(zlist>=model.max_z) = [];
        zlist                           = unique(sort(zlist));
        model                           = fine_model(model, Parameters);

        for k = 1:length(zlist)

            %when bootstrapping, don't use all of them
            if probability < rand

                continue

            end

            mn = model;

            if zlist(k) > upper_z

                mn.(field).z    = [ mn.(field).z; zlist(k)/model.max_z ];
                mn.(field).node = [ mn.(field).node; log(interp1(model.z, model.(field).model, zlist(k), Parameters.interp, exp(model.vs.node(end)))) ];

            else

                mn.(field).z    = [ mn.(field).z; zlist(k)/model.max_z ];
                mn.(field).node = [ mn.(field).node; (log(upper_vs) - 4*((upper_z - zlist(k))/(upper_z/Parameters.division))*Parameters.h) ];%top nodes are made *a little* slower to avoid monotonicity constraint

            end

            [mn.(field).z, ix] = sort(mn.(field).z);
            mn.(field).node    = mn.(field).node(ix);

            mn    = update_z(mn, Parameters);
            mn    = fine_model(mn, Parameters);
            model = mn;

        end

    end

    model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
        
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

