function [ model, f0 ] = fill_J_sed(model, Disp, allWfs, Parameters)

    h = Parameters.h;

    Parameters.surf_96 = false;
    %Parameters.anirec  = true;
    func = @(x) wrapper_res_sed(x, model, Disp, allWfs, Parameters);

    x0 = [ model.vs.node(1); model.vs.node(2); log(model.vs.z(1)*Parameters.plot_z); model.vpvs_block; model.vpvs_z];

    f0 = func(x0);

    J = zeros(length(f0), length(x0));

    sigma = [];

    for k = 1:length(allWfs)

        sigma = [ sigma; allWfs(k).rfr_std/sqrt((allWfs(k).g/Parameters.sample_rate)) ];

    end

    sigma = [ sigma; ones(length(f0) - length(sigma), 1) ];

    if Parameters.anirec

        %first, do it in parallel with the surface waves turned off.
        %surf96 can't be run in parallel
        parfor k = 1:length(x0)
    
            xp1 = x0;
            xn1 = x0;
    
            xp1(k) = xp1(k) + h;
            xn1(k) = xn1(k) - h;
    
            fp1 = func(xp1);
            fn1 = func(xn1);
            J(:, k) = ((0.5*fp1 - 0.5*fn1)/(h));%./sigma;
    
        end

    end

    %model.res = model.res./sigma';

    %then rerun it with surf96 enabled in serial
    %and replace the surface wave gradients
    %get the full thing to get prior information if anirec was skipped
    Parameters.surf_96 = true;
    Parameters.anirec  = false;
    func = @(x) wrapper_res_sed(x, model, Disp, allWfs, Parameters);

    surf_ix = ((length(allWfs)*length(allWfs(1).rfr) + 1):length(model.res))';

    sigma = ones(size(surf_ix));
    sigma(1:length(Disp.c_r), 1) = Disp.c_rstd';

    for k = 1:length(x0)

        xp1 = x0;
        xn1 = x0;

        xp1(k) = xp1(k) + h;
        xn1(k) = xn1(k) - h;

        fp1 = func(xp1);
        fn1 = func(xn1);
        J(surf_ix, k) = ((0.5*fp1(surf_ix) - 0.5*fn1(surf_ix))/(h));%./sigma;

    end

    model.J = J;

end