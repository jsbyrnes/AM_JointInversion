function res = wrapper_res(x, model, Disp, allWfs, Parameters)

    if isrow(x)

        x = x';

    end

    if any(isnan(x))

        res = 1e9*ones(size(model.res'));
        return

    end

    for k = 1:length(Parameters.fields_vec)

        field                  = Parameters.fields_vec{k};
        model.(field).node     = x(1:model.(field).n);
        x(1:(model.(field).n)) = [];

    end

    if model.vpvs_on

        model.vpvs_z     = x(1);
        model.vpvs_block = x(2);
        %model.sed_model  = x(3:7);
        model.max_z      = x(3);

    else

        %model.sed_model  = x(1:5);
        model.max_z      = x(1);

    end

    model       = update_z(model, Parameters);

    model = evaluate_reflectivity(model, allWfs, Disp, Parameters, false)';
    res   = model.res';

end
