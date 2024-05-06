function am = prepare_anirec(model)

    am.vp    = model.vs.model .* model.vpvs.model;
    am.vs    = model.vs.model;
    am.vpvs  = am.vp./am.vs;
    am.z     = model.z;
    am.rho   = nafedrake_rho(am.vp); %g/cm^3
    am.theta = zeros(size(model.z));
    am.phi   = zeros(size(model.z));
    am.A     = zeros(size(model.z));
    am.B     = zeros(size(model.z));
    am.C     = zeros(size(model.z));
    am.t     = model.t;

end