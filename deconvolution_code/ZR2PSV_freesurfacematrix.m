function [ cc ] = ZR2PSV_freesurfacematrix( Z_comp, R_comp, T_comp, slowness, flag, test_ind)
%ZR2PSV [ source_comp_after conv_comp_after ] = ZR2PSV( source_comp_before, conv_comp_before, rp, surfv)
%   Rotate with a free surface transfer matrix. 
%
%   Joseph Byrnes
%   Oct 15th 2012, jbyrnes@uoregon.edu
%   Modified to normalize T on July 21st, 2015.
%

    Vp = 3.50:0.1:7;
    VpVs = 1.5:0.1:2;

    cc = zeros(length(Vp), length(VpVs));
    
    for i = 1:length(Vp)
    
        for j = 1:length(VpVs)
        
            fs_matrix = zeros(3);
            
            qalpha = ( Vp(i)^(-2) - slowness^2)^(0.5);
            qbeta  = ( (Vp(i)/VpVs(j))^(-2) - slowness^2)^(0.5);
            
            fs_matrix(1,1) = ((((Vp(i)/VpVs(j))^2)*slowness^2) - 0.5) / ( Vp(i)*qalpha);
            fs_matrix(1,2) = slowness*((Vp(i)/VpVs(j))^2)/Vp(i);
            fs_matrix(2,1) = slowness*(Vp(i)/VpVs(j));
            fs_matrix(2,2) = (0.5 - (((Vp(i)/VpVs(j))^2)*slowness^2)) / ( (Vp(i)/VpVs(j))*qbeta);
            fs_matrix(3,3) = 0.5;
            
            b = [ Z_comp(:); R_comp(:); T_comp(:) ];
            
            a = fs_matrix*b;
            
            cc(i,j) = corrcoef(a(test_ind,1), a(test_ind,2));
        
        end
   
    end
    
%     [~, index] = min(cc(:));
%         
%     inc_angle = angle(index);
%     
%     R_trace = cosd(inc_angle)*R_comp - sind(inc_angle)*Z_comp;
%     Z_trace = sind(inc_angle)*R_comp + cosd(inc_angle)*Z_comp;
%     
%     if strcmpi(flag, 'P')
%     
%         amp = max(abs(Z_trace));
%         
%     elseif strcmpi(flag, 'SV')
%         
%         amp = max(abs(R_trace));
%         
%     elseif strcmpi(flag, 'SH')
%         
%         amp = max(abs(T_comp));
%         
%     end
%     
%     Z_trace = Z_trace/amp;
%     R_trace = R_trace/amp;
%     T_trace = T_comp/amp;
    
end