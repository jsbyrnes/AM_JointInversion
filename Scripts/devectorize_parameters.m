function model = devectorize_parameters(model, Parameters)

    model.parameterized.sed_z   = model.parameterized.vec(1)*Parameters.MD.sed_z(2)   + Parameters.MD.sed_z(1);
    model.parameterized.moho_z  = model.parameterized.vec(2)*Parameters.MD.moho_z(2)  + Parameters.MD.moho_z(1);
    model.parameterized.moho_dz = model.parameterized.vec(3)*Parameters.MD.moho_dz(2) + Parameters.MD.moho_dz(1); 

    model.parameterized.sed     = model.parameterized.vec(4:4+(length(model.parameterized.sed)-1))...
        *Parameters.MD.sed(2) + Parameters.MD.sed(1);

    model.parameterized.crust   = model.parameterized.vec(4+(length(model.parameterized.sed)):...
        (4+(length(model.parameterized.sed) + length(model.parameterized.crust) - 1)))...
        *Parameters.MD.crust(2) + Parameters.MD.crust(1);

    model.parameterized.mantle  = model.parameterized.vec(((length(model.parameterized.vec) - 3)...
        - length(model.parameterized.mantle) + 1):(length(model.parameterized.vec) - 3))...
        *Parameters.MD.mantle(2) + Parameters.MD.mantle(1);

    model.parameterized.vpvs    = model.parameterized.vec(length(model.parameterized.vec)...
        - 2:end)*Parameters.MD.vpvs(2) + Parameters.MD.vpvs(1);

end