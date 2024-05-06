function surf96model = make_surf96model(model, Parameters)

    %interpolate onto a fixed grid for surf96
    %vpvs = interp1(model.vpvs.z, model.vpvs.node, Parameters.surf_z, 'linear', 'extrap');
    %vs   = interp1(model.vs.z,   model.vs.node,   Parameters.surf_z, 'linear', 'extrap');
    vs   = interp1(model.z,   model.vs.model,   Parameters.surf_z, 'nearest', 'extrap');
    %vs(Parameters.surf_z < model.vs.z(1)) = model.vs.node(1);

%     if any(contains(Parameters.fields_vec, 'vpvs'))
% 
%         vpvs = interp1(model.vpvs.z, model.vpvs.node, Parameters.surf_z, 'linear', model.vpvs.node(end));
%         vpvs(Parameters.surf_z < model.vpvs.z(1)) = model.vpvs.node(1);
% 
%     else
% 
%         vpvs = 1.76*ones(size(vs));
% 
%     end
% 
%     vp  = vs.*vpvs;

    if model.z(1) > 0

        model.z  = [ 0; model.z];
        model.vp = [ model.vp(1); model.vp];

    end

    vp = interp1(model.z, model.vp, Parameters.surf_z, 'nearest', 'extrap');

    rho = nafedrake_rho(vs);%interp1(model.z, model.rho,     Parameters.surf96_z, 'linear', 'extrap');
    %vs  = model.vs.model;

    vp(isnan(vp))   = vp(find(~isnan(vp), 1, 'last'));
    vs(isnan(vs))   = vs(find(~isnan(vs), 1, 'last'));
    rho(isnan(rho)) = rho(find(~isnan(rho), 1, 'last'));

    thick       = diff([ 0; Parameters.surf_z ]);
    %thick        = diff([ 0; model.z ]);

    surf96model = [ thick vp vs rho];
    surf96model = [ surf96model; surf96model(end, :)];

    surf96model(end,1) = 0;

end