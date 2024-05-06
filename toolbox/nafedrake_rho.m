function [ rho ] = nafedrake_rho( vs )
%NAFEDRAKE_RHO Returns the density given by the Nafe-Drake curve. Taken from
% Brocher, 2005 BSSA "Empirical Relations between Elastic Wavespeeds and
% Density in the Crust". Units are vp = [ km/s] and rho = [g/cm^3]. 
% jbyrnes@uoregon.edu, 2012

    if ~ismatrix(vs) || ~isreal(vs)
        
        error('Something wrong with vp - nafedrake_rho.m');
        
    end

    %equation (1)
    %rho = 1.6612*vp - 0.4721*vp.^2 + 0.0671*vp.^3 - 0.0043*vp.^4 + 0.000106*vp.^5;
    rho = 1.227 + 1.53*vs - 0.837*vs.^2 + 0.207*vs.^3 - 0.0166*vs.^4;

end

