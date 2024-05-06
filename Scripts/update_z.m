function model = update_z(model, Parameters)

    %loop through the model and create Parameters.zint per wavelength or
    %1/Parameters.zint as a minimum number of layers
    maxf = max(Parameters.g_array);

%     halfspace_z  = exp(model.sed_model(1));
%     halfspace_vs = exp(model.sed_model(2));
%     exp_decay    = exp(model.sed_model(3));
%     exp_dVs      = exp(model.sed_model(4));
%     gradient     = exp(model.sed_model(5));

    model.vs.z = model.vs.z*model.max_z;

    if isempty(Parameters.zfrac)

        model.z = (0:Parameters.dz_lim(1):Parameters.plot_z)';
        return

    end

    %dz = Parameters.zfrac*exp(model.sed_model(1))/maxf;
    %z  = [ (dz:dz:(halfspace_z))'; halfspace_z ];

    dz = Parameters.zfrac*exp(model.vs.node(1))/maxf;
    z  = [ (dz:dz:(model.vs.z(1)))'; model.vs.z(1) ];

    while z(end) < model.max_z

%         if z(end) < Parameters.sed_model_z
% 
%             if z(end) < halfspace_z
% 
%                 vs = halfspace_vs;
% 
%             else
% 
%                 zp    = z(end) - halfspace_z;
%                 vs = halfspace_vs + gradient*zp + exp_dVs*(1 - 2*(exp(-zp/exp_decay)./(1 + exp(-zp/exp_decay))));
% 
%             end
% 
%         else
% 
%             vs = interp1(model.vs.z, model.vs.node, z(end), Parameters.interp, 'extrap');
% 
%         end

        vs = interp1(model.vs.z, exp(model.vs.node), z(end), Parameters.interp, 'extrap');

        dz = max([ Parameters.dz_lim(1); Parameters.zfrac*vs/maxf ]);
        %z  = [ z; ((model.vs.z(k)+dz):dz:(model.vs.z(k+1)))'; model.vs.z(k+1)];
        z = [ z; z(end) + dz];

    end

    z = sort([ z; model.vs.z]);

%     for k = 1:(length(model.vs.node)-1)
% 
%         dz = max([ Parameters.dz_lim(1); Parameters.zfrac*min([ model.vs.node(k) model.vs.node(k+1) ])/maxf ]);
%         z  = [ z; ((model.vs.z(k)+dz):dz:(model.vs.z(k+1)))'; model.vs.z(k+1)];
% 
%     end

    z = unique(z);

    z(z>=Parameters.plot_z) = [];

    if max(z) < Parameters.plot_z

        model.z = [z; Parameters.plot_z];

    else

        model.z = z;

    end

    model.vs.z = model.vs.z/model.max_z;

end

%     %loop through the model and create Parameters.zint per wavelength or
%     %1/Parameters.zint as a minimum number of layers
%     maxf = max(Parameters.g_array);
%     %wl = interp1(model.vs.z, model.vs.node, 0, 'linear', 'extrap')/maxf;
%     wl = model.vs.node(1)/maxf;
%     z  = wl*((Parameters.zfrac):(Parameters.zfrac):1)';
% 
%     while max(z) < Parameters.plot_z
% 
%         if max(z) > model.vs.z(1)
% 
%             wl = interp1(model.vs.z, model.vs.node, max(z), 'linear', 'extrap')/maxf;
% 
%         else
% 
%             wl = model.vs.node(1)/maxf;
% 
%         end
% 
%         z  = [ z; max(z) + wl*((Parameters.zfrac):(Parameters.zfrac):1)' ];
% 
%     end
