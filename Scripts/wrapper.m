function lpst = wrapper(x, model, Disp, allWfs, Parameters, field,  ind)

    model.(field).node(ind) = x;

    model  = evaluate_reflectivity(model, allWfs, Disp, Parameters, false);
    lpst = model.IC;

end

%attempt to change z at the same time as v 
%     %first, get dz
%     if ind > 1
% 
%         dz = model.(field).z(ind+1) - model.(field).z(ind-1);
% 
%     else
% 
%         dz = model.(field).z(2);
% 
%     end
% 
%     model.(field).z(ind)    = dz*tanh(x(1))/2 + dz/2;%transforms to conitnue
