function model = fill_fields(model)

    if isfield(model, 'vpvs')

        model.vp       = model.vs.model.*(model.vpvs.model);

    else

        model.vp       = model.vs.model.*1.76;

    end

    model.rho      = nafedrake_rho(model.vs.model); %g/cm^3
    
    model.theta = zeros(size(model.z));
    model.phi   = zeros(size(model.z));
    model.A     = zeros(size(model.z));
    model.B     = zeros(size(model.z));
    model.C     = zeros(size(model.z));

end