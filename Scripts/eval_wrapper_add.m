function lpst = eval_wrapper_add(x, model, allWfs, Parameters)

    model.n = model.n + 1;

    model.moveout_t(modeln.n, :)   = x(1:3);
    model.moveout_amp(modeln.n, :) = [ Parameters.prior_amp(1)*randn() Parameters.prior_amp(2)*randn() randn()*Parameters.prior_amp(3) ];
    model.w(modeln.n,1)            = rand()*(log(model.t(end)/2) - 5*log(1/Parameters.sample_rate)) + 5*log(1/Parameters.sample_rate);
    modeln                          = vectorize(modeln);

    model.vec = x;
    model.n   = round(length(model.vec)/7);
    model     = evaluate_optimize(model, allWfs, Parameters);
    lpst      = -1*model.lpst;

end