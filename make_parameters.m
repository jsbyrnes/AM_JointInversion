function Parameters = make_parameters(dirname, fields, values)

    requesting_window = [ -100 250 ];%around the pick
    
    Parameters.save_name = dirname;
    Parameters.tag       = '_060824';

    %for handling the data
    Parameters.sample_rate      = 20;
    Parameters.max_time         = 8;%for the inversion
    Parameters.pre              = 2.5;
    Parameters.datawin          = [ -10 25 ];%for the RFs calculationm
    Parameters.g_array          = 2.^(-2:1);%;2.^(0:3);%logspace(log10(0.5), log10(4), 5);
    Parameters.dt_plot          = 0.05;
    Parameters.tplot            = [-2 25];%wtf is this anymore
    Parameters.downsample       = 1;%fraction of period
    Parameters.p_array          = linspace(0.04, 0.08, 8);
    Parameters.t_array          = [ 5 6 8 10 12 15 20 25 32 40 50 60 ];
    Parameters.limits.vs        = [ 0.5 7];
    Parameters.limits.vpvs      = [ 1.8 0.5 ];%mean/std log
    Parameters.dz_lim           = [ 0.01 10 ];
    Parameters.debug            = false;
    Parameters.plot             = true;
    Parameters.plot_z           = 90;
    Parameters.total_time       = requesting_window(2) - requesting_window(1);
    Parameters.t                = requesting_window(1):(1/Parameters.sample_rate):requesting_window(2);
    Parameters.h                = 0.005;
    Parameters.interp           = 'pchip';
    Parameters.rf_style         = 'multi-taper';
    Parameters.filter_style     = 'butter';
    Parameters.Firstgrad_damp   = [];%log. Positive is more damping.
    Parameters.Secondgrad_damp  = 3;%log. Positive is more damping.
    Parameters.smoothing_style  = 'cauchy';
    Parameters.deconvolve       = true;
    Parameters.IC               = 'BIC';
    Parameters.zfrac            = 1/8;
    Parameters.division         = 3;%for making new nodes. 
    Parameters.zstd             = 5;%for scaling the depth of nodes
    Parameters.sed_model_z      = 6;%in km. When inversions are allowed
    Parameters.baz_range        = [];
    Parameters.monotonicity_w   = 5;
    Parameters.min_error        = 3;%in percent of peak amplitude
    Parameters.bootstrapping    = false;
    Parameters.bts_samples      = 50;
    Parameters.hd_error         = 0.33;
    Parameters.baz              = [ ];

    Parameters.move_out.perform_moveout = true;
    Parameters.move_out.ref_p           = 0.06;
    %I made this correction by hand for crustal structure
    Parameters.move_out.amp_slope       = 18;%constant correction
    Parameters.move_out.time_slope      = 0.85;%correction per s of lag

    Parameters.requesting_window = requesting_window;%around the pick
    Parameters.surf_z   = (linspace(0.2, sqrt(Parameters.plot_z), 100)').^2;
    Parameters.fields_vec = { 'vs' };%Parameters.fields_vec = { 'vs', 'vpvs'};
    
    Parameters.vpvs_block = true;

    if nargin==3

        for k = 1:length(fields)

            Parameters.(fields{k}) = values{k};

        end

    end

    if Parameters.move_out.perform_moveout

        Parameters.p_array = Parameters.move_out.ref_p;

    end

    Parameters.anirec  = true;
    Parameters.surf_96 = true;

end