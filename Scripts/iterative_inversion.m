function [ model, allWfs ] = iterative_inversion(Parameters, rawData, Disp, truemodel)

    if ~isempty(Parameters.baz)

        if Parameters.baz(2) > Parameters.baz(1)

            ind = [rawData(:).baz]>Parameters.baz(1) & [rawData(:).baz] < Parameters.baz(2);

        else

            ind = [rawData(:).baz]>Parameters.baz(2) & [rawData(:).baz] < Parameters.baz(1);

        end

    else

        ind = 1:length(rawData);

    end

    allWfs       = fill_RFs(2.5, 2.5*1.76, Parameters.t, rawData(ind), Parameters, 500);

    model = node_inversion(Parameters, allWfs, Disp, truemodel, []);
    dv    = 100;
    
    while abs(dv) > 0.05
    
        vs0            = model.vs.model(1);
        disp([ '*** Reinverting with surface velocity of ' num2str(vs0) ])
    
        allWfs         = fill_RFs(vs0, vs0*exp(model.vpvs_block), Parameters.t, rawData, Parameters, 500);
        model          = node_inversion(Parameters, allWfs, Disp, truemodel, model);
    
        dv = model.vs.model(1) - vs0;
        disp([ ' -> Surface velocity changed by ' num2str(dv) ])

    end

end

%     vs_vec = 0.25:0.25:2.5;
%     
%     vs0            = vs_vec(1);
%     disp([ 'Inverting with surface velocity of ' num2str(vs0) ])
% 
%     allWfs         = fill_RFs(vs0, vs0*1.76, Parameters.t, rawData, Parameters, 1);
%     model          = node_inversion(Parameters, allWfs, Disp, truemodel, [], vs0);
% 
%     for k = 2:length(vs_vec)
% 
%         vs0          = vs_vec(k);        
%         disp([ 'Inverting with surface velocity of ' num2str(vs0) ])
%     
%         allWfs         = fill_RFs(vs0, vs0*1.76, Parameters.t, rawData, Parameters, 1);
%         model(k)       = node_inversion(Parameters, allWfs, Disp, truemodel, model(1), vs0);
%     
%     end
% 
%     keyboard





%     g0      = Parameters.g_array;
%     model   = [];
% 
%     for k = 1:length(Parameters.g_array)
%     
%         Parameters.g_array = g0(1:k);
% 
%         if isempty(model)
% 
%             allWfs       = fill_RFs(1, 1*1.76, Parameters.t, rawData, Parameters, 500);
% 
%         else
% 
%             allWfs       = fill_RFs(model.vs.model(1), model.vs.model(1)*exp(model.vpvs_block), Parameters.t, rawData, Parameters, 500);
% 
%         end
% 
%         disp([ '********** First inversion with frequencies up to ' num2str(Parameters.g_array(end)) ' Hz' ])
%         dv = 100;
%         
%         model    = node_inversion(Parameters, allWfs, Disp, truemodel, model);
%         
%         while abs(dv) > 0.05
%         
%             vs0            = model.vs.model(1);
%             disp([ '*** Reinverting with surface velocity of ' num2str(vs0) ])
%         
%             allWfs         = fill_RFs(vs0, vs0*exp(model.vpvs_block), Parameters.t, rawData, Parameters, 500);
%             model          = node_inversion(Parameters, allWfs, Disp, truemodel, model);
%         
%             dv = model.vs.model(1) - vs0;
%             disp([ ' -> Surface velocity changed by ' num2str(dv) ])
% 
%         end
%     
%     end
