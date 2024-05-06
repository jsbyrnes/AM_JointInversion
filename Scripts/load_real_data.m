function [ rawData, Disp, Parameters ] = load_real_data(dirname, Parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Real Data! Body wave part
        
    files = dir('./data/*.itr');
    
    x = strsplit(dirname, '.');
    
    network = x{1};
    staname = x{2};
        
    sta = 4;
    lta = 20;
        
    total = 0;

    %used later
    sta_lat = [];
    sta_lon = [];

    disp([ '*******Starting download for ' staname])

    for k = 1:length(files)
    
        disp(['Making RF #' num2str(k) ' of ' num2str(length(files)) ])
    
        [D,T0,H] = rdsac([ './data/' files(k).name ]);
        
        phase_list = 'P';
        %get the predicted travel time
        urlstr     = sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&noheader=true&mintimeonly=true&phases=%s',H.GCARC,H.EVDP/1000,phase_list);
            
        count = 0;

        while count < 25

            try
            
                tStr  = urlread(urlstr);
                tCell = textscan(tStr,'%f %f %s %f %f %f %f %f %s %s');
            
                count = 26;

            catch
        
                pause(5)
                count = count + 1;
    
                if count == 25

                    disp('-> Error getting phase times')

                end
    
                continue
        
            end
    
        end

        %phases= tCell{3};
        times = tCell{4};
    
        %make the time of the RF to get
        ds = strsplit(H.KZDATE);
    
        ds = [ ds{1} '-' ds{2} '-' ds{4} ' ' H.KZTIME];
    
        originTimeNum  = datenum(ds, 'mm-dd-yyyy HH:MM:SS.FFF');
        
        startTime = originTimeNum + (times + Parameters.requesting_window(1))/(24*60*60); 
        endTime   = originTimeNum + (times + Parameters.requesting_window(2))/(24*60*60);
    
        startTime = datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');   
        endTime   = datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');
    
        myTrace = irisFetch.Traces(network,...
            staname,'*','BHZ,BHN,BHE',startTime, endTime);
    
        if length(unique({ myTrace(:).location })) ~= 1

            loc_vec  = { myTrace(:).location };
            chan_vec = { myTrace(:).channel };

            locations = unique(loc_vec);

            for lix = 1:length(locations)

                %check if all three components are there
                z_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BHZ'));
                n_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BHN'));
                e_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BHE'));

                if any(z_check) && any(n_check) && any(e_check)

                    myTrace = [ myTrace(z_check) myTrace(n_check) myTrace(e_check) ];
                    break

                end

            end

        end

        if length(myTrace)~=3 || length(unique({ myTrace(:).location })) ~= 1
    
            %try again with BH1 and BH2
            myTrace = irisFetch.Traces(network,...
                staname,'*','BHZ,BH1,BH2',startTime, endTime);
        
            if length(unique({ myTrace(:).location })) ~= 1
    
                loc_vec  = { myTrace(:).location };
                chan_vec = { myTrace(:).channel };
    
                locations = unique(loc_vec);
    
                for lix = 1:length(locations)
    
                    %check if all three components are there
                    z_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BHZ'));
                    n_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BH1'));
                    e_check = (strcmp(loc_vec, locations{lix}) & strcmp(chan_vec, 'BH2'));
    
                    if any(z_check) && any(n_check) && any(e_check)
    
                        myTrace = [ myTrace(z_check) myTrace(n_check) myTrace(e_check) ];
                        break
    
                    end
    
                end
    
            end

            if length(myTrace)~=3
            
                disp(' -> Some traces not loaded')
                continue

            else

                %make them north/east
                names = {myTrace.channel};

                %find the north channel (might be prerotation)
                oneind = find(contains(names, 'BH1'));
            
                %find the east channel (might be prerotation)
                twoind = find(contains(names, 'BH2'));

                d1 = myTrace(oneind).data;
                d2 = myTrace(twoind).data;

                try

                    myTrace(oneind).data  = cosd(myTrace(oneind).azimuth)*d1 - sind(myTrace(oneind).azimuth)*d2;
                    myTrace(twoind).data  = sind(myTrace(oneind).azimuth)*d1 + cosd(myTrace(oneind).azimuth)*d2;

                catch

                    %sample count might be uneven
                    disp(' -> Some traces not loaded')
                    continue

                end

                myTrace(oneind).channel = 'BHN'; 
                myTrace(twoind).channel = 'BHE'; 

            end
            
        end

        if length(unique([ myTrace(:).sampleCount ])) ~= 1

            disp(' -> Sample count uneven')
            continue
            
        end

        names = {myTrace.channel};
    
        % Get rid of empty structures or missing data
        data  = {myTrace.data};
        tf_empty = cellfun('isempty',data);
        myTrace  = myTrace(~tf_empty);
    
        %find the north channel (might be prerotation)
        vind = find(contains(names, 'BHZ'));
        
        %find the north channel (might be prerotation)
        nind = find(contains(names, 'BHN'));
    
        %find the east channel (might be prerotation)
        eind = find(contains(names, 'BHE'));
    
        if length(vind) > 1
    
            [~, ind] = max([myTrace(nind).sampleRate]);
            vind = vind(ind);
    
            [~, ind] = max([myTrace(nind).sampleRate]);
            nind = nind(ind);
    
            [~, ind] = max([myTrace(eind).sampleRate]);
            eind = eind(ind);
    
        end
        
        if isempty(sta_lon)

            sta_lon = myTrace(1).longitude;
            sta_lat = myTrace(1).latitude;
            
        end

        slow(k) = H.USER2;
        baz(k)  = H.BAZ;

        myTrace = wfResample_jsb(myTrace, Parameters);    

        myTrace(vind).data = bandpassfilt_rfs(myTrace(vind).data, 1/Parameters.sample_rate, 9.5, 1/50);
        myTrace(nind).data = bandpassfilt_rfs(myTrace(nind).data, 1/Parameters.sample_rate, 9.5, 1/50);
        myTrace(eind).data = bandpassfilt_rfs(myTrace(eind).data, 1/Parameters.sample_rate, 9.5, 1/50);

        R       = myTrace(nind).data*cosd(H.BAZ - 180)  + myTrace(eind).data*sind(H.BAZ - 180);
        T       = -myTrace(nind).data*sind(H.BAZ - 180) + myTrace(eind).data*sind(H.BAZ - 180);

        myTrace(nind).data = R;
        myTrace(eind).data = T;

        Z = myTrace(vind).data;%bandpassfilt_rfs(myTrace(vind).data, 1/Parameters.sample_rate, 9.5, 1/50);
        %Z = bandpassfilt_rfs(myTrace(vind).data, 1/Parameters.sample_rate, 9.5, 1/50);
        %R = bandpassfilt_rfs(myTrace(nind).data, 1/Parameters.sample_rate, 9.5, 1/50);
        
        %This happens for bad traces
        if rms(R)/rms(Z) > 3 || max(abs(R))/max(abs(Z)) > 2

            continue

        end

        Zpicking = bandpassfilt_rfs(myTrace(vind).data, 1/Parameters.sample_rate, 1.5, 0.5);
        
%         if strcmp(Parameters.filter_style, 'butter')
% 
%             Z = bandpassfilt_rfs(myTrace(vind).data, 1/Parameters.sample_rate, 5, 1/50);
%             R = bandpassfilt_rfs(myTrace(nind).data, 1/Parameters.sample_rate, 5, 1/50);
% 
%         elseif strcmp(Parameters.filter_style, 'gaussian')
% 
%             gsig = 1/(pi*7.5);
%             gauss_sig = round(gsig/0.05);
%     
%             x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
%             Gauss_win = exp(-x.^2/(2*gauss_sig^2));
%             Z = conv(myTrace(vind).data,Gauss_win,'same');
%             R = conv(myTrace(nind).data,Gauss_win,'same');
% 
%         end

        tpicking = (0:(length(Z)-1))/Parameters.sample_rate + Parameters.requesting_window(1);
    
        [pick, snr] = sta_lta(Zpicking, 0, 33, ...
            sta, lta, tpicking);
    
        [~, ind] = min(abs(tpicking - (pick - sta/2)));
    
        disp(snr)

        %only keep 150 s, wont all be used
        try
    
            Z = Z((ind + Parameters.datawin(1)*Parameters.sample_rate):(ind + 150*Parameters.sample_rate));
            R = R((ind + Parameters.datawin(1)*Parameters.sample_rate):(ind + 150*Parameters.sample_rate));
    
        catch
    
            disp(' -> Trace length incorrect. Pick probably too early')
            continue%very rare
    
        end
    
        if snr < 30
    
            disp(' -> Signal to noise ratio too low')            
            continue

        end

        total = total + 1;

        if isrow(Z)

            Z = Z';

        end

        if isrow(R)

            R = R';

        end

        taper = (abs(Parameters.datawin(1))/(0.5*(length(Z)/Parameters.sample_rate)));

        rawData(total).Z     = Z.*tukeywin(length(Z), taper);
        rawData(total).R     = R.*tukeywin(length(Z), taper);
        rawData(total).p     = slow(k);
        rawData(total).baz   = baz(k);
        rawData(total).rfr   = [];
        rawData(total).rqual = 1;%doesnt really help
                    
    end

    disp([ num2str(total) ' traces loaded out of ' num2str(length(files)) ' possible' ])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %surface wave phase velocities (now that we know the station location!)
    %lots of stuff is hard wired!
    for k = 1:length(Parameters.t_array)

        if Parameters.t_array(k) < 20

            %load from Ekstrom
            data = readmatrix([ '../obs_dispersion/R' num2str(Parameters.t_array(k)) ...
                '_USANT15.pix'], 'FileType', 'text', 'NumHeaderLines', 11);

            %really should be read from the file headers but this was
            %easier
            switch Parameters.t_array(k)
            
                case 5

                    ref = 2.9833;

                case 6

                    ref = 3.0532;

                case 8

                    ref = 3.1467;

                case 10

                    ref = 3.2140;

                case 12

                    ref = 3.2760;

                case 15

                    ref = 3.3677;

            end

            Disp.c_r(k) = ref*(1 + 0.01*griddata(data(:, 1), data(:, 2), data(:, 4), sta_lat, sta_lon));

        else

            %load from Jin
            data = readmatrix([ '../obs_dispersion/helmholtz_stack_LHZ_' num2str(Parameters.t_array(k)) ...
                '.xyz'], 'FileType', 'text');

            Disp.c_r(k) = griddata(data(:, 1), data(:, 2), data(:, 3), sta_lat, sta_lon);
            
        end

    end
    
    Disp.c_r = Disp.c_r';
    Disp.c_rstd = abs(140 ./ (140 / 4 + Parameters.t_array' / 100) - 4)/2;
    Disp.c_rstd(Disp.c_rstd < 0.025) = 0.025;

    nans = isnan(Disp.c_r);

    if ~isempty(nan)

        Parameters.t_array(nans) = [];
        Disp.c_r(nans) = [];
        Disp.c_rstd(nans) = [];

    end

end
