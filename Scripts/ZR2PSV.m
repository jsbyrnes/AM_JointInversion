function [ P, SV ] = ZR2PSV( Z_comp, R_comp, alpha, beta, p)
%Kennet, 1991 equation 18

    P  = ((p*beta^2)/alpha)*R_comp + (((beta^2)*p^2 - 0.5)/(alpha*sqrt(1/alpha^2 - p^2)))*Z_comp;
    SV = ((0.5 - (beta^2)*p^2)/(beta*sqrt(1/beta^2 - p^2)))*R_comp + p*beta*Z_comp;

end

%     angle = -80:1:80;
% 
%     for i = 1:length(angle)
%     
%         R_trace_tmp = cosd(angle(i))*R_comp - sind(angle(i))*Z_comp;
%         %Z_trace_tmp = sind(angle(i))*R_comp + cosd(angle(i))*Z_comp;
%                     
%         test(i) = rms(R_trace_tmp(window(1):window(2)));
%    
%     end
%     
%     [~, index] = min(test);
%         
%     R_trace = cosd(angle(index))*R_comp - sind(angle(index))*Z_comp;
%     Z_trace = sind(angle(index))*R_comp + cosd(angle(index))*Z_comp;
%     
%     if strcmpi(flag, 'P')
%     
%         amp = max(abs(Z_trace));
%         
%     else
%         
%         amp = max(abs(R_trace));
%         
%     end
%     
%     Z_trace = Z_trace/amp;
%     R_trace = R_trace/amp;
