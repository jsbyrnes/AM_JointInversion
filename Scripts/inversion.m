function model = inversion(model, Parameters, allWfs, Disp, truemodel, maxiter)
    
    if nargin == 5

        maxiter = 100;

    end

    opt_single = optimset('MaxIter', 1e9, 'MaxFunEvals', 1e9, 'TolX', 0.01);   
    func       = @(x) wrapper_res(x, model, Disp, allWfs, Parameters);

    x0 = [];

    for k = 1:length(Parameters.fields_vec)

        x0 = [ x0; model.(Parameters.fields_vec{k}).node ];

    end

    if model.vpvs_on

        x0 = [ x0; model.vpvs_z; model.vpvs_block ];

    end
    %x0 = [ x0; model.sed_model];
    x0 = [ x0; model.max_z ]; 

    nvar = length(x0);

    dx = 1;
    df = 1;

    neval = 0;
    f0    = func(x0);
    disp([ ' -> Starting inversion at a ' Parameters.IC ' of ' num2str(model.IC) ]);

    count = 0;

    lambda  = 0.01;
    l_stall = 0;

    while count < maxiter && df > 0.25
        
        count = count + 1;

        ic = model.IC;

        x0 = [];
        for k = 1:length(Parameters.fields_vec)
    
            x0 = [ x0; model.(Parameters.fields_vec{k}).node ];
    
        end

        if model.vpvs_on

            x0 = [ x0; model.vpvs_z; model.vpvs_block ];

        end
        %x0 = [ x0; model.sed_model];
        x0 = [ x0; model.max_z ];

        model = fill_J(model, Disp, allWfs, Parameters);
        
        %f2 = @(lambda) sum( func(x0 - inv(model.J'*model.J + exp(lambda)*diag(diag(model.J'*model.J)))*model.J'*f0).^2 );
        %f2 = @(lambda) sum( func(x0 - inv(model.J'*model.J + exp(lambda)*eye(length(x0)))*model.J'*f0).^2 );

        %[lambda, ~, ~, ~] = fminbnd(f2, -5, 10, opt_single);

        fn = 1e9;

        while sum(fn.^2) > sum(f0.^2) && lambda < 1e3

            if l_stall > 3

                lambda = lambda/2;%maybe you can speed it back up?
                l_stall = 0;
                
            end

            dx    = inv(model.J'*model.J + lambda*diag(diag(model.J'*model.J)))*model.J'*f0;
            %dx    = inv(model.J'*model.J + exp(lambda)*eye(length(x0)))*model.J'*f0;
            %dx = sign(dx).*min(abs(dx), 0.5*ones(size(dx)));
            
            if (x0(nvar-1) - dx(nvar-1)) < log(1.1)%min Vp/Vs considered

                dx(nvar-1) = x0(nvar-1) - log(1.1);

            end
            
            fn = func(x0 - dx);

            if sum(fn.^2) > sum(f0.^2)

                lambda  = lambda*2;
                l_stall = 0;

            else

                l_stall = l_stall + 1;

            end

        end

        if sum(fn.^2) < sum(f0.^2) && lambda < 1e3

            f0 = fn;
            x0 = x0 - dx;

            for k = 1:length(Parameters.fields_vec)
    
                field                        = Parameters.fields_vec{k};
                model.(field).node           = x0(1:model.(field).n);
                x0(1:(model.(field).n))    = [];

            end

            if model.vpvs_on
        
                model.vpvs_z     = x0(1);
                model.vpvs_block = x0(2);

                %model.sed_model  = x0(3:7);

                model.max_z      = x0(3);

            else

                %model.sed_model  = x0(1:5);
                model.max_z      = x0(1);
                
            end

            model       = update_z(model, Parameters);
            model       = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);

            df = model.IC - ic;

            plot_model_new(5678, model, Parameters, allWfs, Disp, truemodel);
            %plot_data(101, Parameters, model.synWfs, model.Disp);

            if isempty(lambda)
    
                disp([ '   -> ' Parameters.IC ' increased to ' num2str(model.IC) ]);
    
            else
    
                disp([ '   -> ' Parameters.IC ' increased to ' num2str(model.IC) '. λ at ' num2str(lambda) ]);
    
            end

        else

            df = 0;

            disp('   -> Search failed to improve')

        end

    end
    
    model = update_z(model, Parameters);
    model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
    
end

%attempted hierarchical inversion, but not needed!
%     %update the errors
%     func = @(x) wrapper_sig(x, model, Disp, allWfs, Parameters);
%     opt  = optimset('MaxIter', 1e100, 'MaxFunEvals', 1e100);
%     %vec  = fminsearch(func, [model.sig_cr model.sig_cl model.sig_rf], opt);%NM handles bounds well
%     vec  = fminsearch(func, -3*[1 1 ones(size(model.sig_rf))], opt);%NM handles bounds well
% 
%     model.sig_cr = vec(1);
%     model.sig_cl = vec(2);
%     model.sig_rf = vec(3:end);
%     model        = evaluate_reflectivity(model, allWfs, Disp, Parameters, true);
% 
%     model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false); %update the data

%     %%%%%%now relocate on the L curve
%     neval = 0;
%     dx    = 1;
%     h     = 0.05;
% 
%     disp(' -> Adjusting the smoothing parameter')
% 
%     while max(abs(dx)) > 1e-2
% 
%         model_fine = fine_model(model, Parameters);
%     
%         func = @(x) wrapper_res(x, model, Disp, allWfs, Parameters, false);%without smoothing
%         x0 = [ model.vs.node; model.vpvs.node ];
%         %f0 = func(x0);
% 
%         %first, get the Jacobin at current values
%         J = zeros(length(allWfs(1).rf)*numel(allWfs) + length(Disp.c_r) + length(model_fine.z), length(x0));
%         for k = 1:length(x0)
%     
%             xp = x0;
%             xn = x0;
%     
%             xp(k) = xp(k) + h;
%             xn(k) = xn(k) - h;
%     
%             J(:, k) = (func(xp) - func(xn))/(2*h);
%             neval = neval + 2;
% 
%         end
% 
%         %newton's search for best smoothing parameters a la best lpst
%         dk    = 1;
%         count = 0;
% 
%         while abs(dk) > 0.1 && count < 5
% 
%             smth = model.curvature + (-5:0.25:5);
%     
%             for k = 1:length(smth)
%     
%                 Jn = [ J; smoothing_J(model, Parameters, smth(k), h);];
%     
%                 mn           = model;
%                 mn.curvature = smth(k);
%                 %get f0 at this smoothing to invert against
%                 func_smth = @(x) wrapper_res(x, mn, Disp, allWfs, Parameters, true);%with new smoothing
%                 f0smth    = func_smth(x0);
%     
%                 %now, get the fit at the current values
%                 %f2 = @(lambda) sum( func_smth(x0 - inv(Jn'*Jn + exp(lambda)*eye(length(x0)))*Jn'*f0smth).^2 );
%                 %[lambda, fval, ~, out] = fminbnd(f2, -3, log(1000), opt_single);
%         
%                 %dx       = inv(Jn'*Jn + exp(lambda)*eye(length(x0)))*Jn'*f0smth;%get dx considering smoothing
%                 %big change here - this is just the Gauss-Newton step. LM
%                 %step might be better, but then you need to evaluate the
%                 %dampening term and that's hard to do! but can be done if
%                 %needed, just commented out
% 
%                 dx       = (Jn'*Jn)\Jn'*f0smth;%get dx considering smoothing
%                 %fsmth(k) = sum(func_smth(x0 - dx).^2);%get lpst
%                 dx_smth(k) = log(sum(dx.^2));
%                 fsmth(k) = log(sum(func_smth(x0 - dx).^2));%just the fit (with priors on normal paramters)
% 
%                 neval = neval + 2;% + out.funcCount;
%     
%             end
% 
%             if unique(sign(diff(fsmth)))==-1%check for monotonicity
% 
%                 nf     = min(fsmth):0.01:max(fsmth);
%                 ndx    = interp1(fsmth, dx_smth, nf, 'spline');
% 
%                 %3rd and 4th derivative on fit not counting smoothing. This is the
%                 %first and second derivative of the "curvature" ie the second
%                 %derivative, which I want to maximize
%                 %coefficients from https://en.wikipedia.org/wiki/Finite_difference_coefficient
%                 %d3 = (-0.5*fsmth(1) + fsmth(2) - fsmth(4) + 0.5*fsmth(5))/0.25^3;
%                 %d4 = (fsmth(1) - 4*fsmth(2) + 6*fsmth(3) - 4*fsmth(4) + fsmth(5))/0.25^4;
%     
%                 [~, ix] = min(abs(nf - fsmth(3)));
% 
%                 d3 = (-0.5*ndx(ix-2) + ndx(ix-1) - ndx(ix+1) + 0.5*ndx(ix+2))/0.01^3;
%                 d4 = (ndx(ix-2) - 4*ndx(ix-1) + 6*ndx(ix) - 4*ndx(ix+1) + ndx(ix+2))/0.25^4;
% 
%                 if sign(d4) == 1 %wrong curvative, just do gradient ascent
%     
%                     df = 0.25*sign(d3);
%     
%                 else %should be good
%     
%                     df = -d3/d4;% in the direction of the gradient when curvative is negative
%     
%                     %d = sign(dk)*min([ 0.75 abs(dk) ]);
%     
%                 end
%     
%                 step_f = fsmth(3) + df;
% 
%                 %and what smooth is this?
%                 model.curvature = interp1(fsmth, smth, step_f, 'spline');
%                 
%                 %model.curvature = model.curvature + dk;
% 
%             else
% 
%                 nsmth = min(smth):0.01:max(smth);
%                 nf    = interp1(smth, fsmth, nsmth, 'spline');
%                 [~, ix] = min(nf);
% 
%                 model.curvature = nsmth(ix);
%                 break
% 
%             end
% 
%             count = count + 1;
% 
%         end
% 
%         Jn = [ J; smoothing_J(model, Parameters, model.curvature, h);];
%     
%         func_smth = @(x) wrapper_res(x, model, Disp, allWfs, Parameters, true);%with new smoothing
%         f0smth    = func_smth(x0);
% 
%         dx = (Jn'*Jn)\Jn'*f0smth;%get dx considering smoothing
% 
%         %and get dx at this smoothing
%         fn = func(x0 - dx);%note that the fit might get worse! That's ok
%         f0 = fn;
%         x0 = x0 - dx;
% 
%         disp([ '   -> Fit now ' num2str(sum(f0.^2)) ' with max step of ' ...
%             num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls.' ]);
% 
% %             if isempty(lambda)
% %     
% %                 disp([ '   -> Fit reduced to ' num2str(sum(f0.^2)) ' with max step of ' ...
% %                     num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls.' ]);
% %     
% %             else
% %     
% %                 disp([ '   -> Fit reduced to ' num2str(sum(f0.^2)) ' and max step of ' ...
% %                     num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls. λ raised to ' num2str(exp(lambda)) ]);
% %     
% %             end
% 
% %         else
% % 
% %             disp('   -> Fit cannot improve without over fitting')
% % 
% %         end
% 
%         model.vs.node   = x0(1:model.vs.n);
%         model.vpvs.node = x0(model.vs.n+1 : end);
%         model           = evaluate_reflectivity(model, allWfs, Disp, Parameters, false); %update the data
% 
%     end
% 
% 
% 
% function Jsmth = smoothing_J(model, Parameters, smth, h)
% 
%     jix  = 1;
% 
%     model = fine_model(model, Parameters);
% 
%     %get the derivatives with the smoothing
%     for fidx = 1:length(Parameters.fields_vec)
% 
%         field = Parameters.fields_vec{fidx};
%         
%         for j = 1:model.(field).n
% 
%             resp = [];
%             resn = [];
% 
%             mnp = model;
%             mnn = model;
%             mnp.(field).node(j) = mnp.(field).node(j) + h;
%             mnn.(field).node(j) = mnn.(field).node(j) - h;
% 
%             mnp      = fine_model(mnp, Parameters);
%             mnn      = fine_model(mnn, Parameters);
% 
%             dz = mnp.z(2) - mnp.z(1);
%             
%             %need to change limits to be defined by a structure
%             %with matching field names
% 
%             if strcmp(field, 'vs')%non-dimensionalization different for vs (uniform distribution) than for others 
% 
% %                 resp = [ gradient(gradient((mnp.(field).model' - Parameters.limits.(field)(1))/diff(Parameters.limits.(field)), dz), dz)/exp(smth) ...
% %                     zeros(1, (length(Parameters.fields_vec)-1)*length(model.z)) ];        
% %                 resn = [ gradient(gradient((mnn.(field).model' - Parameters.limits.(field)(1))/diff(Parameters.limits.(field)), dz), dz)/exp(smth) ...
% %                     zeros(1, (length(Parameters.fields_vec)-1)*length(model.z))];        
%                 resp = [ gradient((mnp.(field).model' - Parameters.limits.(field)(1))/diff(Parameters.limits.(field)), dz)/exp(smth) ...
%                     zeros(1, (length(Parameters.fields_vec)-1)*length(model.z)) ];        
%                 resn = [ gradient((mnn.(field).model' - Parameters.limits.(field)(1))/diff(Parameters.limits.(field)), dz)/exp(smth) ...
%                     zeros(1, (length(Parameters.fields_vec)-1)*length(model.z))];        
% 
%             else
% 
%                 resp = [ zeros(1, (fidx-1)*length(model.z)) ...
%                     gradient((mnp.(field).model' - Parameters.limits.(field)(1))/Parameters.limits.(field)(2), dz)/exp(smth) ...
%                      zeros(1, (length(Parameters.fields_vec) - fidx)*length(model.z))];        
%                 resn = [ zeros(1, (fidx-1)*length(model.z)) ...
%                     gradient((mnn.(field).model' - Parameters.limits.(field)(1))/Parameters.limits.(field)(2), dz)/exp(smth) ...
%                      zeros(1, (length(Parameters.fields_vec) - fidx)*length(model.z))];        
%                 
%             end
% 
%             Jsmth(:, jix) = ((resp - resn)/(2*h))';
%             jix = jix + 1;
% 
%         end
% 
%     end
% 
% end

%     %%%%%%now relocate on the L curve
%     neval = 0;
%     dx    = 1;
%     
%     disp(' -> Adjusting the smoothing parameter')
% 
%     while max(abs(dx)) > 1e-2
% 
%         model_fine = fine_model(model, Parameters);
%     
%         func = @(x) wrapper_res(x, model, Disp, allWfs, Parameters, false);%without smoothing
%         x0 = [ model.vs.node; model.vpvs.node ];
%         f0 = func(x0);
% 
%         %first, get the Jacobin at current values
%         J = zeros(length(allWfs(1).rf)*numel(allWfs) + length(Disp.c_r) + length(model_fine.z), length(x0));
%         for k = 1:length(x0)
%     
%             xp = x0;
%             xn = x0;
%     
%             xp(k) = xp(k) + h;
%             xn(k) = xn(k) - h;
%     
%             J(:, k) = (func(xp) - func(xn))/(2*h);
%             neval = neval + 2;
% 
%         end
% 
%         %newton's search for best smoothing parameters
%         dk    = 1;
%         count = 0;
% 
%         while abs(dk) > 0.1 && count < 5
% 
%             smth = model.curvature + (-0.5:0.25:0.5);
%     
%             for k = 1:length(smth)
%     
%                 Jn = [ J; smoothing_J(model, Parameters, smth(k), h);];
%     
%                 mn           = model;
%                 mn.curvature = smth(k);
%                 %get f0 at this smoothing to invert against
%                 func_smth = @(x) wrapper_res(x, mn, Disp, allWfs, Parameters, true);%with new smoothing
%                 f0smth    = func_smth(x0);
%     
%                 %now, get the fit at the current values
%                 %f2 = @(lambda) sum( func_smth(x0 - inv(Jn'*Jn + exp(lambda)*eye(length(x0)))*Jn'*f0smth).^2 );
%                 %[lambda, fval, ~, out] = fminbnd(f2, -3, log(1000), opt_single);
%         
%                 %dx       = inv(Jn'*Jn + exp(lambda)*eye(length(x0)))*Jn'*f0smth;%get dx considering smoothing
%                 %big change here - this is just the Gauss-Newton step. LM
%                 %step might be better, but then you need to evaluate the
%                 %dampening term and that's hard to do! but can be done if
%                 %needed, just commented out
% 
%                 dx       = (Jn'*Jn)\Jn'*f0smth;%get dx considering smoothing
%                 fsmth(k) = sum(func(x0 - dx).^2);%now evaluate without considering smoothing
%         
%                 neval = neval + 2;% + out.funcCount;
%     
%             end
% 
%             %3rd and 4th derivative on fit not counting smoothing. This is the
%             %first and second derivative of the "curvature" ie the second
%             %derivative, which I want to maximize
%             %coefficients from https://en.wikipedia.org/wiki/Finite_difference_coefficient
%             d3 = (-0.5*fsmth(1) + fsmth(2) - fsmth(4) + 0.5*fsmth(5))/0.25^3;
%             d4 = (fsmth(1) - 4*fsmth(2) + 6*fsmth(3) - 4*fsmth(4) + fsmth(5))/0.25^4;
% 
%             if sign(d4) == 1 %wrong curvative, just do gradient ascent
% 
%                 dk = 0.25*sign(d3);
% 
%             else %should be good
% 
%                 dk = -d3/d4;% in the direction of the gradient when curvative is negative
% 
%                 dk = sign(dk)*min([ 0.75 abs(dk) ]);
% 
%             end
% 
%             model.curvature = model.curvature + dk;
% 
%             count = count + 1;
% 
%         end
% 
%         Jn = [ J; smoothing_J(model, Parameters, model.curvature, h);];
%     
%         func_smth = @(x) wrapper_res(x, model, Disp, allWfs, Parameters, true);%with new smoothing
%         f0smth    = func_smth(x0);
% 
%         dx = (Jn'*Jn)\Jn'*f0smth;%get dx considering smoothing
% 
%         %and get dx at this smoothing
%         fn = func(x0 - dx);%note that the fit might get worse! That's ok
%         f0 = fn;
%         x0 = x0 - dx;
% 
%         disp([ '   -> Fit now ' num2str(sum(f0.^2)) ' with max step of ' ...
%             num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls.' ]);
% 
% %             if isempty(lambda)
% %     
% %                 disp([ '   -> Fit reduced to ' num2str(sum(f0.^2)) ' with max step of ' ...
% %                     num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls.' ]);
% %     
% %             else
% %     
% %                 disp([ '   -> Fit reduced to ' num2str(sum(f0.^2)) ' and max step of ' ...
% %                     num2str(max(abs(dx))) ' with ' num2str(neval) ' function calls. λ raised to ' num2str(exp(lambda)) ]);
% %     
% %             end
% 
% %         else
% % 
% %             disp('   -> Fit cannot improve without over fitting')
% % 
% %         end
% 
%         model.vs.node   = x0(1:model.vs.n);
%         model.vpvs.node = x0(model.vs.n+1 : end);
%         model           = evaluate_reflectivity(model, allWfs, Disp, Parameters, false); %update the data
% 
%     end







%         dx    = inv(J'*J)*J'*f0;
%         fn    = func(x0 - dx);
%         neval = neval + 1;
% 
%         if (sum(f0.^2) - sum(fn.^2)) < 0.1 || max(abs(dx)) < 0.01
% 
%             f2 = @(lambda) sum( func(x0 - inv(J'*J + exp(lambda)*eye(length(x0)))*J'*f0).^2 );
% 
%             [lambda, ~, ~, out] = fminbnd(f2, -3, log(1000), opt_single);
% 
%             dx    = inv(J'*J + exp(lambda)*eye(length(x0)))*J'*f0;
%             fn    = func(x0 - dx);
% 
%             neval = neval + 1 + out.funcCount;
% 
%         end

%         while sum(fn.^2) > sum(f0.^2)
% 
%             if lambda == 0
%                 
%                 lambda = 0.01;
% 
%             else
% 
%                 lambda = lambda*5;
% 
%             end
% 
%             dx    = inv(J'*J + lambda*eye(length(x0)))*J'*f0;
%             fn    = func(x0 - dx);
% 
%         end

%     opt      = optimoptions('lsqnonlin', 'FiniteDifferenceStepSize', 5e-2, 'FiniteDifferenceType', 'central', ...
%         'FunctionTolerance', 0.001, 'Display', 'iter', 'StepTolerance', 1e-3, ...
%         'Algorithm', 'levenberg-marquardt', 'InitDamping', 1, 'ScaleProblem', 'jacobian');
% 
%     func = @(x) wrapper_res(x, model, Disp, allWfs, Parameters);
%     vec  = lsqnonlin(func, [ model.vs.node; model.vpvs.node ], [ Parameters.vs_lim(1)*ones(size(model.vs.node)); 0.5*ones(size(model.vpvs.node)) ], ...
%         [ Parameters.vs_lim(2)*ones(size(model.vs.node)); 5*ones(size(model.vpvs.node)) ], opt);

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

% old code to do 1 d optimization in parallel. Works... ok, kinda shitty
%         lp = log10(1000);
%         ln = log10(0.001);
% 
%         while lp - ln > 0.1
% 
%             lsearch = logspace(ln, lp, Parameters.numworkers);
% 
%             parfor k = 1:length(lsearch)
% 
%                 fsearch(k) = f2(lsearch(k));
% 
%             end
% 
%             neval = neval + Parameters.numworkers;
% 
%             [~, ind] = min(fsearch);
% 
%             if ind == 1
% 
%                 dl = log10(lsearch(2)) - log10(lsearch(1));
% 
%                 ln = log10(lsearch(1)) - 2*dl;
%                 lp = log10(lsearch(1)) + 2*dl;
% 
%             elseif ind == Parameters.numworkers || length(unique(fsearch))==1
% 
%                 dl = log10(lsearch(Parameters.numworkers)) - log10(lsearch(7));
% 
%                 ln = log10(lsearch(Parameters.numworkers)) - 2*dl;
%                 lp = log10(lsearch(Parameters.numworkers)) + 2*dl;
% 
%             else 
%             
%                 ln = log10(lsearch(ind-1));
%                 lp = log10(lsearch(ind+1));
% 
%             end
% 
%         end
%         lv      = linspace(lsearch(1), lsearch(end), 1e4);
%         [~, ix] = min(interp1(lsearch(fsearch<1e18), fsearch(fsearch<1e18), lv, 'spline'));
%         lambda  = lv(ix);
%         %lambda = 10^((lp + ln)/2);

