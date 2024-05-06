function lpst = wrapper_sig(x, model, Disp, allWfs, Parameters, ix)

    %model.sig_cr = x(1);
    %model.sig_cl = x(2);
    %model.sig_rf = x(3:(3+length(Parameters.g_array)-1));
    model.sig_rf(ix) = x;

    model = evaluate_reflectivity(model, allWfs, Disp, Parameters, true);
    lpst  = -1*model.lpst;

end
