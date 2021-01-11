function par=set_parameters_master(par)

% TEMPLATE MATCHING
par.max_spk = 20000;                % max. # of spikes before starting templ. match.

% FORCE MEMBERSHIP PARAMETERS
par.force_auto = 'y';               %automatically force membership
par.template_sdnum = 2.5;            % max radius of cluster in std devs.
par.template_k = 10;                % # of nearest neighbors
par.template_k_min = 10;            % min # of nn for vote
%par.template_type = 'mahal';        % nn, center, ml, mahal
par.template_type = 'center';        % nn, center, ml, mahal
par.force_feature = 'spk';          % feature use for forcing (whole spike shape)
%par.force_feature = 'wav';         % feature use for forcing (wavelet coefficients).

% FEATURES PARAMETRES
%par.features = 'wav';              %choice of spike features
%par.features = 'wavfpv'            % wavelets features and a forced peak-valley value are used for sorting
par.features = 'all'                % wavelets, PCA and peak-valley value are used for sorting


par.feature_points=48;              % Number of points used to calculate the features, this used to be 24
par.inputs = 10;                    %number of inputs to the clustering
par.scales = 4;                     %scales for wavelet decomposition
if strcmp(par.features,'pca');      %number of inputs to the clustering for pca
    par.inputs=3; 
end

%SPC PARAMETERS
par.mintemp = 0;                    %minimum temperature
par.maxtemp = 0.201;                %maximum temperature
par.tempstep = 0.01;                %temperature step
par.num_temp = floor(...
(par.maxtemp - ...
par.mintemp)/par.tempstep);         % total number of temperatures 
par.stab = 0.8;                     % stability condition for selecting the temperature
par.SWCycles = 100;                 % number of montecarlo iterations
par.KNearNeighb = 11;               % number of nearest neighbors
par.randomseed = 0;                 % if 0, random seed is taken as the clock value
%par.randomseed = 147;              % If not 0, random seed   
par.fname_in = 'tmp_data';          % temporary filename used as input for SPC
par.fname = 'data';                 % filename for interaction with SPC
par.max_spikes = par.max_spk/2;       % maximum number of spikes to plot
par.min_clus_abs = 20;              %minimum cluster size (absolute value)
par.min_clus_rel = 0.001;            %minimum cluster size (relative to the total nr. of spikes)
par.temp_plot = 'log';              %temperature plot in log scale
par.max_clus = 7;                   % maximum number of clusters, must be at least 
par.max_histos=8;


%GET SPIKES prameters
par.w_pre=par.feature_points/2;     %number of pre-event data points stored
par.w_post=par.feature_points/2;    %number of post-event data points stored
par.interpolation = 'y';            %interpolation for alignment
par.int_factor = 2;                 %interpolation factor
par.detect_fmin = 300;              %high pass filter for detection
par.detect_fmax = 3000;             %low pass filter for detection
par.sort_fmin = 300;                %high pass filter for sorting
par.sort_fmax = 3000;               %low pass filter for sorting
par.segments = 1;                   %nr. of segments in which the data is cutted.
par.tmax='all';

% HISTOGRAM PARAMETERS
for i=1:par.max_histos+1
    eval(['par.nbins' num2str(i-1) ' = 10;']);  % # of msec for the range of the ISI histograms
    eval(['par.bin_step' num2str(i-1) ' = 0.2;']);  % bin step in msec
    eval(['par.violation_mark' num2str(i-1) ' = 1;']);  % msec where the violation line is drawn
end

%Correct multichannel artifact
par.correct_multichannel_artifact=1;
par.exc_dt=0.001;

%For behavior
par.doBehavior=1; %lick bits are not read in if this is set to 0
par.dstep=500;    %Data are resampled to reduce the sampling rate by this number
par.no_below=6;   %Time of decision is assigned to the first point below p=0.05 that 
                    %has no_below trailing points <0.05
                    
par.read_entire_file=1;




