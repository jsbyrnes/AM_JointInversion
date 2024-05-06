function model = devectorizeMF(model)

    vec = model.vec;

    model.moveout_t    = reshape(vec(1:(model.n*6)), [ model.n 6 ]);
    vec(1:(model.n*6)) = [];
    model.moveout_amp  = reshape(vec(1:(model.n*6)), [ model.n 6 ]);
    vec(1:(model.n*6)) = [];
    model.w  = reshape(vec(1:(model.n*6)), [ model.n 6 ]);
    vec(1:(model.n*6)) = [];
    model.sig          = vec(end);

end