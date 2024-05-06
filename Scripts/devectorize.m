function model = devectorize(model)

    vec = model.vec;

    model.moveout_t    = reshape(vec(1:(model.n*3)), [ model.n 3 ]);
    vec(1:(model.n*3)) = [];
    model.moveout_amp  = reshape(vec(1:(model.n*3)), [ model.n 3 ]);
    vec(1:(model.n*3)) = [];
    model.w            = vec(1:(model.n));
    model.sig          = vec(end);

end