function lpst = eval_wrapper(x, model, allWfs, Parameters)

    model.vec = x;
    model.n   = round((length(model.vec)-1)/7);
    model     = evaluate_optimize(model, allWfs, Parameters);
    lpst      = -1*model.lpst;

end