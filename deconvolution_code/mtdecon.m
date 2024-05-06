function [ output_args ] = mtdecon( C )
%MTDECON do the multitaper deconvolution
% SMH 10/10

% general parameters
%==========================================================================

rf_norm  = 1;          % normalize the RFs to the maximum in a window around the direct arrival? (1 = yes, 0 = no)
norm_win = [-3 3];     % window in sec around the direct arrival

%==========================================================================

%get stations
istn =  [C.p2data 'META.station.' C.db_name '.mat'];

load(istn);

% record control parameters
RfC.decon_style     = C.decon_style;
RfC.rf_filter      = [];
RfC.filter_nc      = [];
RfC.rf_shift       = C.rf_shift;
RfC.rf_norm        = rf_norm;
RfC.norm_win       = norm_win;

% print out parameters
fprintf('\nmake_rfs.m\n==========\n')

fprintf('\nDB name     : %s',C.db_name)
fprintf('\nInput path  : %s',C.p2postcull)
fprintf('\nOutput path : %s\n',C.p2rf)

fprintf('\nDecon method: %s\n',C.decon_style)

if use_source_est
    fprintf('\nUsing source estimate for decon.')
    fprintf('\nSpectral smoothing: %s',C.smooth_style)
    fprintf('\nSpectral estimator: %s\n',C.spec_est_style)
else
    fprintf('\nUsing direct component trace for decon.')
end

if C.rf_shift < 1.5*C.taper_time
    fprintf('\n\n***RF shift time less than 1.5x taper time!\nIncreasing shift time...\n')
    C.rf_shift     = C.taper_time*1.5;
    RfC.rf_shift = C.rf_shift;
end
fprintf('\nRF time shift: %2.2f sec\n',C.rf_shift)

fprintf('\nNormalize RFs     : %d',rf_norm)
if rf_norm
    fprintf('\nNormalizing window: %2.3f to %2.3f sec about direct arrival\n',norm_win(1),norm_win(2))
end

% check paths
if ~exist(C.p2postcull,'dir')
    error('path: %s is not a valid directory!',C.p2postcull)
end

if ~exist(C.p2rf,'dir')
    s = unix(['mkdir ' C.p2rf]);
    if s
        error('unable to make RF dir!')
    end
end

% DECON TRACES, go by station
%==========================================================================

[nsta, ~] = size(S);

% load SOURCE file
if use_source_est
    sfile = [ C.p2source 'SOURCE.' C.db_name '.' C.smooth_style '.' C.spec_est_style  '.mat' ];
    load(sfile)
else
    sfile = [];
end

for ii=1:nsta   % loop over stations in B struct
    
    fprintf('\nDecon: station %d of %d\n',ii,nsta)
    fprintf('=============================\n')
    
    nbin = B(ii).nbin;
    sid  = B(ii).sid;
    
    fprintf('Station sid: %d\n',sid)
    fprintf('\n# of bins: %d\n',nbin)
    
    % init RF struct
    RF_D = struct('C1',[],'C2',[],'C3',[]);
    % init
    T_master = [];
    S0 = [];
    
    RMS  = zeros(3,nbin);
    MAMP = zeros(3,nbin);
    FR   = zeros(3,nbin);
    
    phase_type = cell(1,nbin);
    
    for jj=1:nbin % loop over event bins for this structure
        
        eid_list   = B(ii).eid_list{jj};        % list of events in bin
        phase_list = B(ii).phase_list{jj};      % list of phases in bin
        
        nevt = length(eid_list);
        
        % subset Source structure and check the event distribution in B vs
        % Source (if requested)
        if use_source_est
            % init
            % ipx=[];
            ip = zeros(nevt,1);
            for kk=1:nevt % loop over events in bin
                % find event in source struct
                e_ip = find([Source.eid]==eid_list(kk));
                % check eid
                if isempty(e_ip) % incongruency between source_spec_est and binning
                    fprintf('\n\n***event eid: %d not found in source structure!\n',eid_list(kk))
                    keyboard
                    % fprintf('Skipping\n')
                    % ipx=[ipx kk]; % ignore this event
                end
                
                % find phase
                p_ip=false(length(e_ip),1);
                for ll=1:length(e_ip)
                    temp_ip  = strmatch(phase_list{kk},Source(e_ip).phase,'exact');
                    p_ip(ll) = ~isempty(temp_ip);
                end
                
                if ~any(p_ip)
                    fprintf('\n\n***Phase: %s not found in source structure for event eid: %d !\n',phase_list{kk},eid_list(kk))
                    keyboard
                elseif sum(p_ip)~=1
                    fprintf('\n\n***Phase: %s found in source structure multiple times for event eid: %d !\n',phase_list{kk},eid_list(kk))
                    keyboard
                else
                    ip(kk) = e_ip(p_ip); % record
                end
                
            end % event loop
            
            Source_cut = Source(ip); % subset/arrange source
        end
        
        % grab the data for this station-bin
        [T D] = collect_trace_data3(C,sid,eid_list,phase_list);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
        
        % make some checks
        if T.n3Ctra ~= nevt
            fprintf('\n\n***Number of traces requested different than read!\n')
            keyboard
        end
        
        if ~all(T.eid == eid_list')   % shouldn't happen
            fprintf('\n\n***Read eid list is different than requested!\n')
            keyboard
        end
        
        if length(char(T.phase)) ~= length(char(phase_list)) % shouldn't happen
            fprintf('\n\n***Read phase list is different than requested!\n')
            keyboard
        else
            if ~all(char(T.phase) == char(phase_list)) % shouldn't happen
                fprintf('\n\n***Read phase list is different than requested!\n')
                keyboard
            end
        end
        
        % check some more things (consistency in the structures)
        npts_s = unique([Source_cut.npts]);
        npts_t = unique([T.npts]);
        
        sps_t = unique([T.sps]);
        sps_s = unique([Source_cut.sps]);
        
        df = unique([Source.df]);
        
        if length(npts_s) ~= 1 || length(npts_t) ~= 1 || length(sps_t) ~= 1 || length(sps_s) ~= 1 || length(df) ~= 1
            fprintf('\n\n***Change in the time parameters detected in T and/or Source structures, sid: %d, bid: %d \n',sid,B(ii).bid(jj))
            keyboard
        end
        
        if npts_s ~= npts_t && ~C.zpad
            fprintf('\n\n***Differing number of time samples in T and Source structures (1)!\n')
            keyboard
        elseif npts_s ~= 2*npts_t && C.zpad
            fprintf('\n\n***Differing number of time samples in T and Source structures (2)!\n')
            keyboard
        end
                
        if sps_s ~= sps_t
            fprintf('\n\n***Differing number of samples per second in T and Source structures!\n')
            keyboard
        end
        
        npts = npts_t;
        sps  = sps_t;
        
%         [b,a] = butter(filter_nc,rf_filter.*(2*1/sps)); % filter coeffs

        % look up component ips
        [cmp_ip rot_style] = look_up_cmp_ip([T.cmp_name]);
        rot_style = unique(rot_style);
        
        if length(rot_style) > 1
            fprintf('\n\n***Mixed trace rotation detected, sid: %d, bid: %d !\n',sid,B(ii).bid(jj))
            keyboard
        end
        
        % check the phase
        lst_u = last_upper([T.phase]);
        lst_u = unique(lst_u);
                
        phase_type(jj) = lst_u(1); % keep track of phase type for this bin
        
        if length(lst_u) > 1
            fprintf('\n\n***Mixed arrival phase detected, sid: %d, bid: %d !\n',sid,B(ii).bid(jj))
            keyboard
        end
        
        % arrange trace data
        switch lst_u{1}
            case 'S'
                [C2 C1 C3] = collect_trace_cmp(D,cmp_ip); % C1=SV,C2=P,C3=SH
                % fix cmp_ip to reflect this switching
                cmp_ip = [cmp_ip(2,:);cmp_ip(1,:);cmp_ip(3,:)];
            case 'P'
                [C1 C2 C3] = collect_trace_cmp(D,cmp_ip); % C1=P,C2=SV,C3=SH
            otherwise
                fprintf('\n\n***Something weird with the phase type for event %d!\n',eid)
                keyboard
        end
        
        % check traces
        if ~any(any(C1,1)) || ~any(any(C2,1)) %|| ~any(any(C3,1))
            fprintf('\n\n***zero trace detected!\n')
            keyboard
        end
        
        if ~all(all(isfinite(C1))) || ~all(all(isfinite(C2))) || ~all(all(isfinite(C3)))
            fprintf('\n\n***NaN/Inf trace value detected!\n')
            keyboard
        end
        
        % record arranging
        T.cmp_name = T.cmp_name(cmp_ip);
        
        if C.zpad
            z = zeros([npts nevt]);
            C1 = [C1;z];
            C2 = [C2;z];
            C3 = [C3;z];
            
            T.npts = 2*npts; % need to remember the zero padding in the T struct
            npts = 2*npts;
        end
        
        switch C.decon_style
            case 'MPLSQ'
                
                % parameter
                %==============================================================
                filt_norm = 1;    % use allpass or pseudo allpass, i.e. normalize the phase filters? (yes = 1, no = 0, .5 = semi-normalize)
                %==============================================================
                
                % write control flag
                RfC.filt_norm = filt_norm;
                
                if use_source_est
                    % LSQR decon
                    [G So] = LSQdecon2(C1,C2,C3,Source_cut,T,RfC);
                else
                    % LSQR decon
                    [G So] = LSQdecon2(C1,C2,C3,[],T,RfC);
                end
                
                S0 = [S0 ; {So}];
                
            %case 'ETMT'
                
            %case 'ITDD'
                
            %case 'SPARSE'
                                
                
            otherwise
                fprintf('\n\n***Unknown string: %s, in variable''C.decon_style''\n',C.decon_style)
                keyboard
        end
        
        
        
        %******************************************************************
        
        % cut down RF
        if C.zpad
            
            npts = npts/2;
            
            G.C1 = detrend(G.C1(1:npts,1));
            G.C2 = detrend(G.C2(1:npts,1));
            G.C3 = detrend(G.C3(1:npts,1));
            
            T.npts=npts; % need to remember this in the T struct
            
        end
        
        %******************************************************************
        
        % RF taper
        r = C.taper_time/(1/sps*.5*(npts));
        tap = tukeywin(npts,r);
        
        T_master = [T_master ; T];
        
%         % filter
%         G.C1 = single(filtfilt(b,a,double(G.C1.*tap).*tap)); % hardwired single here*************
%         G.C2 = single(filtfilt(b,a,double(G.C2.*tap).*tap));
%         G.C3 = single(filtfilt(b,a,double(G.C3.*tap).*tap));
        
        G.C1 = G.C1.*tap;
        G.C2 = G.C2.*tap;
        G.C3 = G.C3.*tap;
        
        % normalize
        if rf_norm
            
            nip = (C.rf_shift+norm_win)*sps;
            nip(nip<1)    = 1;
            nip(nip>npts) = npts;
            
            norm = max(abs(G.C1(nip(1):nip(2))));
            
            G.C1 = G.C1./norm;
            G.C2 = G.C2./norm;
            G.C3 = G.C3./norm;
        end
        
        % calculate quality criteria
        RMS(1,jj) = sqrt(mean(G.C1.^2));
        RMS(2,jj) = sqrt(mean(G.C2.^2));
        RMS(3,jj) = sqrt(mean(G.C3.^2));
        
        MAMP(1,jj) = max(abs(G.C1));
        MAMP(2,jj) = max(abs(G.C2));
        MAMP(3,jj) = max(abs(G.C3));
        
        df = 1/(npts*1/sps);
        f_ip = round(C.rf_filter./df); % use the filter window for freq win
        
        FR(1,jj) = max(abs(G.FFT_C1(f_ip(1):f_ip(2)))) / median(abs(G.FFT_C1(f_ip(1):f_ip(2)))); % unfiltered data
        FR(2,jj) = max(abs(G.FFT_C2(f_ip(1):f_ip(2)))) / median(abs(G.FFT_C2(f_ip(1):f_ip(2))));
        FR(3,jj) = max(abs(G.FFT_C3(f_ip(1):f_ip(2)))) / median(abs(G.FFT_C3(f_ip(1):f_ip(2))));
        
        % record RFs
        RF_D.C1 = [RF_D.C1 G.C1];
        RF_D.C2 = [RF_D.C2 G.C2];
        RF_D.C3 = [RF_D.C3 G.C3];
        
    end % bin loop
    
    % % remove empty bins from B
    % B(ii).bid(ipbx)=[];
    % B(ii).center_x(ipbx)=[];
    % B(ii).center_y(ipbx)=[];
    % B(ii).mean_delt(ipbx)=[];
    % B(ii).mean_baz(ipbx)=[];
    % B(ii).mean_rayp(ipbx)=[];
    % B(ii).ntra(ipbx)=[];
    
    % B(ii).eid_list(ipbx)=[];
    % B(ii).phase_list(ipbx)=[];
    
    B(ii).nbin = length(B(ii).bid);
    
    % check
    if B(ii).nbin~=size(RF_D.C1,2)
        fprintf('\n\n***# of bins different than # of RFs\n')
        keyboard
    end
    
    % save RF data
    RF_B = B(ii);
    RF_T = T_master;
    
    % create hybrid meta data struct for RF data
    RF_T = BT_to_RFT3(RF_B,RF_T,phase_type,sfile,RfC);
    
    % record quality criteria
    RF_T.RMS = RMS;
    RF_T.max_amp = MAMP;
    RF_T.freq_rat = FR;
    
    RF_T.IBAD = false(3,nbin); % init a culling flag for later
    
    phase_type_str = unique(phase_type);
    
    if length(phase_type_str) ~= 1
        fprintf('\n\n***More than on phase type detected!\n')
        keyboard
    end
    
    rf_file = ['RF.' C.db_name '.' char(phase_type_str) '.' num2str(sid) '.' C.decon_style '.mat'];
    
    fprintf('\nSaving RF data file: %s\n',rf_file)
    
    save([C.p2rf rf_file],'RF_T','RF_D','S0','RfC','SrcC');
    
end % station loop

fprintf('\nOutput path : %s\n',C.p2rf)

fprintf('\nDone.\n')


end

