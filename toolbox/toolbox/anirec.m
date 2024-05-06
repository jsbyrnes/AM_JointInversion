function [ z_comp, r_comp, t_comp ] = anirec( phase, dt, slow, baz, model, rotate )
%ANIREC Runs the anirec program, and returns a normalized three component green's
%function in ray coordinates. Length of returned trace is fixed to 100 seconds. Phase is case
%insensitive. dt in seconds, slow in sec/km, baz in degree clockwise from
%north. Both velocities in km/s, z is depth to base of layer. Last layers
%is infinite half space. 

    phase = upper(phase);

    if ~strcmpi(phase, 'P') && ~strcmpi(phase, 'SV') && ~strcmpi(phase, 'SH')
        
        error('Phase should be P, SV, or SH');
        
    end
        
    model.vp = model.vp*1000;
    model.vs = model.vp./model.vpvs;
    model.z = model.z*1000;
    model.rho = model.rho*1000;

    cc = 1/slow;
    
    if strcmpi(phase, 'P')        

        %called in this fashion so that there is no output to the screen
        evalc('[z_comp, r_comp, t_comp] = anirec_gateway(model.theta, model.phi, model.z, model.vp, model.A, model.B, model.vs, model.C, model.rho, length(model.vp) - 1, cc, baz, dt, 1)');
                
        if rotate
            
            [~, ind] = max(z_comp);         

            try

                [z_comp, r_comp] = ZR2PSV_emp(z_comp, r_comp, (ind - 2):(ind+2));

            catch

                keyboard

            end
            %[z_comp, t_comp] = ZR2PSV(z_comp, t_comp, 'P');
            
        end
        
    elseif strcmpi(phase, 'SV')
        
        [z_comp, r_comp, t_comp] = anirec_gateway(model.theta, model.phi, model.z, model.vp, ...
            model.A, model.B, model.vs, model.C, model.rho, length(model.vp) - 1, cc, baz, dt, 2);
        
        if max(abs(r_comp)) > max(r_comp)
            
            r_comp = r_comp*-1;
            z_comp = z_comp*-1;
            
        end
        
        %z_comp = -1*z_comp;

        if rotate
        
            %[z_comp, r_comp] = ZR2PSV(z_comp, r_comp, 'S');
            [z_comp, r_comp] = ZR2PSV_surfv(z_comp, r_comp, t_comp, 'S', model.vs(1)/1000, slow);
            
        end
        
    elseif strcmpi(phase, 'SH')
                
        [z_comp, r_comp, t_comp] = anirec_gateway(model.theta, model.phi, model.z, model.vp, ...
            model.A, model.B, model.vs, model.C, model.rho, length(model.vp) - 1, cc, baz, dt, 3);
        
        if max(abs(r_comp)) > max(r_comp)
            
            r_comp = r_comp*-1;
            z_comp = z_comp*-1;
            
        end
        
        if rotate
        
            [z_comp, r_comp, t_comp] = ZR2PSV(z_comp, r_comp, t_comp, 'SH');
            
        end
                
    end
    
end

