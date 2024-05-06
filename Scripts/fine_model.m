function model = fine_model(model, Parameters)

    z = model.z;

    %%%%%%%%%%
    %build the vs model
    vs = zeros(size(z));
      
    model.vs.z = model.vs.z*model.max_z;

%     %dimensionalize sediment model
%     halfspace_z  = exp(model.sed_model(1));
%     halfspace_vs = exp(model.sed_model(2));
%     exp_decay    = exp(model.sed_model(3));
%     exp_dVs      = exp(model.sed_model(4));
%     gradient     = exp(model.sed_model(5));
% 
%     model.vs.z = [ Parameters.sed_model_z; model.vs.z ];

%     zp    = Parameters.sed_model_z - halfspace_z;
%     model.vs.node = [ halfspace_vs + gradient*zp + exp_dVs*(1 - 2*(exp(-zp/exp_decay)./(1 + exp(-zp/exp_decay)))); ...
%         model.vs.node ];

    for k = 1:length(vs)
            
        if z(k) < model.vs.z(1)

            vs(k) = exp(model.vs.node(1));

        else

            vs(k) = interp1(model.vs.z, exp(model.vs.node), z(k), Parameters.interp, exp(model.vs.node(end)));

        end
    
%         if z(k) < Parameters.sed_model_z
% 
%             if z(k) < halfspace_z
% 
%                 vs(k) = halfspace_vs;
% 
%             else
% 
%                 zp    = z(k) - halfspace_z;
%                 vs(k) = halfspace_vs + gradient*zp + exp_dVs*(1 - 2*(exp(-zp/exp_decay)./(1 + exp(-zp/exp_decay))));
% 
%             end
% 
%         else
% 
%             vs(k) = interp1(model.vs.z, model.vs.node, z(k), Parameters.interp, 'extrap');
% 
%         end

        %vs(k) = interp1(model.vs.z, model.vs.node, z(k), Parameters.interp, 'extrap');

    end
    model.vs.z     = model.vs.z/model.max_z;
    model.vs.model = vs;

    %remove the fake node
%     model.vs.z    = model.vs.z(2:end);
%     model.vs.node = model.vs.node(2:end);
    %%%%%%%%

    if any(contains(Parameters.fields_vec, 'vpvs'))

        %%%%%%%%
        %build the vpvs model
        vpvs = zeros(size(z));
        model.vpvs.z = model.vpvs.z*model.max_z;
        for k = 1:length(vpvs)
        
            %vpvs(k) = interp1(model.vpvs.z, model.vpvs.node, z(k), 'linear', 'extrap');
    
            if z(k) < model.vpvs.z(1)
    
                vpvs(k) = model.vpvs.node(1);
    
            else
    
                vpvs(k) = interp1(model.vpvs.z, model.vpvs.node, z(k), Parameters.interp, model.vs.node(end));
    
            end
        
        end

        model.vpvs.z = model.vpvs.z/model.max_z;
        model.vpvs.model = vpvs;
        %%%%%%%%%%

    end
    
    model.z = z;

%     %cut down
%     current_ind = 1;
% 
%     while current_ind < (length(model.vs.model)-1)
% 
%         dvs   = model.vs.model(current_ind + 1)   - model.vs.model(current_ind);
%         %dvpvs = model.vpvs.model(current_ind + 1) - model.vpvs.model(current_ind);
% 
%         if abs(dvs) < 2*eps% && abs(dvpvs) < 2*eps
% 
%             model.vs.model(current_ind)   = [];
%             %model.vpvs.model(current_ind) = [];
%             model.z(current_ind)          = [];
% 
%         else
% 
%             current_ind = current_ind + 1;
% 
%         end
% 
%     end

    model   = fill_fields(model);
    
    if model.vpvs_on

        vpvs_z    = exp(model.vpvs_z);
        model.vp  = model.vs.model.*(1.76 + (exp(model.vpvs_block) - 1.76)*exp(-model.z/vpvs_z));

    end

end
