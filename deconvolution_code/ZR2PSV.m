function [ Z_trace, R_trace, T_trace, inc_angle ] = ZR2PSV( Z_comp, R_comp, T_comp, flag)
%ZR2PSV [ source_comp_after conv_comp_after ] = ZR2PSV( source_comp_before, conv_comp_before, rp, surfv)
%   Rotate the Z and R traces into ray coordinates, then properly normalize T.
%
%   Joseph Byrnes
%   Oct 15th 2012, jbyrnes@uoregon.edu
%   Modified to normalize T on July 21st, 2015.
%

    angle = 0:1:90;

    for i = 1:length(angle)
    
        R_trace_tmp = cosd(angle(i))*R_comp - sind(angle(i))*Z_comp;
        Z_trace_tmp = sind(angle(i))*R_comp + cosd(angle(i))*Z_comp;
        
        if strcmpi(flag, 'P')
            
            [~, index] = max(abs(Z_trace_tmp));
            
            test(i) = sum(abs(R_trace_tmp(index - 30:index + 30)));
            
        elseif strcmpi(flag, 'SV') || strcmpi(flag, 'SH')
            
            [~, index] = max(abs(R_trace_tmp));
            
            test(i) = sum(abs(Z_trace_tmp(index - 30:index + 30)));
                            
        end
   
    end
    
    [~, index] = min(test);
        
    inc_angle = angle(index);
    
    R_trace = cosd(inc_angle)*R_comp - sind(inc_angle)*Z_comp;
    Z_trace = sind(inc_angle)*R_comp + cosd(inc_angle)*Z_comp;
    
    if strcmpi(flag, 'P')
    
        amp = max(abs(Z_trace));
        
    elseif strcmpi(flag, 'SV')
        
        amp = max(abs(R_trace));
        
    elseif strcmpi(flag, 'SH')
        
        amp = max(abs(T_comp));
        
    end
    
    Z_trace = Z_trace/amp;
    R_trace = R_trace/amp;
    T_trace = T_comp/amp;
    
end
