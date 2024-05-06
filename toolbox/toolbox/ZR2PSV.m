function [ Z_trace, R_trace ] = ZR2PSV( Z_comp, R_comp, flag)
%ZR2PSV [ source_comp_after conv_comp_after ] = ZR2PSV( source_comp_before, conv_comp_before, rp, surfv)
%   Rotate the Z and R traces into ray coordinates
%   surfv should be the vp(1) for a P and vs(1) for an S
%   Make sure that rp is plane wave(s/km). Source component is the 
%   channel that the incoming wave in on, conv component is the channe;
%   that the converted phases are on. Also renormalizes after rotation
%   to incoming component. Assumes all later phases are less than incoming
%   pulse(probably true for synthetic data, want to be true for real data).
%   Flag should be 'P' or 'S', depending on the incoming wave. Only used
%   for normalization.
%
%   Joseph Byrnes
%   Oct 15th 2012, jbyrnes@uoregon.edu
%

    angle = -80:1:80;

    for i = 1:length(angle)
    
        R_trace_tmp = cosd(angle(i))*R_comp - sind(angle(i))*Z_comp;
        Z_trace_tmp = sind(angle(i))*R_comp + cosd(angle(i))*Z_comp;
        
        if strcmpi(flag, 'P')
            
            [~, index] = max(abs(Z_trace_tmp));
            
            try
            
                test(i) = sum(abs(R_trace_tmp(index - 30:index + 30)));
                
            catch
                
                keyboard
                
            end
            
        else
            
            [~, index] = max(abs(R_trace_tmp));
            
            try
            
                test(i) = sum(abs(Z_trace_tmp(index - 30:index + 30)));
                
            catch
                
                keyboard
                
            end
            
        end
   
    end
    
    [~, index] = min(test);
        
    R_trace = cosd(angle(index))*R_comp - sind(angle(index))*Z_comp;
    Z_trace = sind(angle(index))*R_comp + cosd(angle(index))*Z_comp;
    
    if strcmpi(flag, 'P')
    
        amp = max(abs(Z_trace));
        
    else
        
        amp = max(abs(R_trace));
        
    end
    
    Z_trace = Z_trace/amp;
    R_trace = R_trace/amp;
     
end