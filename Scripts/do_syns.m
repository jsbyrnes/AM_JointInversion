function synWfs = do_syns(model, Parameters)
    
    modeltmp      = model;
    modeltmp.vs   = model.vs.model;

    if contains(Parameters.fields_vec, 'vpvs')

        modeltmp.vpvs = model.vpvs.model;

    else
        
        modeltmp.vpvs = model.vp./model.vs.model;

    end

    parfor k = 1:length(Parameters.p_array)

        slow = Parameters.p_array(k);

        %evalc('[P_comp, SV_comp, ~] = anirec("P", 0.05, slow, 0, model, 0)');
        [Z_comp, R_comp, T_comp] = anirec("P", 0.05, slow, 0, modeltmp, 0);

        %Z_comp = Z_comp - median(Z_comp);
        %R_comp = R_comp - median(R_comp);
        %T_comp = T_comp - median(T_comp);
        
        [~, ind] = max(Z_comp);

        indn   = round(Parameters.pre/0.05) + 1;
        Z_comp = circshift(Z_comp, indn - ind);
        R_comp = circshift(R_comp, indn - ind);
        T_comp = circshift(T_comp, indn - ind);

        f = @(tau) centering(Z_comp, 1/0.05, tau, indn);

        tau = fminbnd(f, -0.05, 0.05);

        Z_comp = delay_continuous(Z_comp, 1/0.05, tau);
        R_comp = delay_continuous(R_comp, 1/0.05, tau);
        T_comp = delay_continuous(T_comp, 1/0.05, tau);
    
%         Z_comp(1:round(5/0.05)) = 0;
%         R_comp(1:round(5/0.05)) = 0;
%         T_comp(1:round(5/0.05)) = 0;
%         Z_comp(round(25/0.05):end) = 0;
%         R_comp(round(25/0.05):end) = 0;
%         T_comp(round(25/0.05):end) = 0;

        synWfs(k).Z = [ Z_comp; zeros(size(Z_comp)) ];
        synWfs(k).R = [ R_comp; zeros(size(R_comp)) ];
        synWfs(k).T = [ T_comp; zeros(size(T_comp)) ];

%         synWfs(k).Z = [ zeros(size(Z_comp)); Z_comp; zeros(size(Z_comp)) ];
%         synWfs(k).R = [ zeros(size(R_comp)); R_comp; zeros(size(R_comp)) ];
%         synWfs(k).T = [ zeros(size(T_comp)); T_comp; zeros(size(T_comp)) ];
%         synWfs(k).Z = Z_comp;
%         synWfs(k).R = R_comp;
%         synWfs(k).T = T_comp;
        synWfs(k).p     = slow;
        synWfs(k).baz   = rand*360;
        synWfs(k).rqual = 10;

    end

end