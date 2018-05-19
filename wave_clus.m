function varargout = wave_clus(varargin)
% WAVE_CLUS M-file for wave_clus.fig
%      WAVE_CLUS, by itself, creates a new WAVE_CLUS or raises the existing
%      singleton*.
%
%      H = WAVE_CLUS returns the handle to a new WAVE_CLUS or the handle to
%      the existing singleton*.
%
%      WAVE_CLUS('Property','Value',...) creates a new WAVE_CLUS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to wave_clus_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      WAVE_CLUS('CALLBACK') and WAVE_CLUS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in WAVE_CLUS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wave_clus


% This is a modification of the wave_clus code from Quian Quiroga
% for use in the Restrepo lab
% 
% The Quian Quiroga code is available in
%
% http://www2.le.ac.uk/centres/csn/research-2/spike-sorting
%
% We have added the ablity to use peak to valley and PCA in addition to the 
% wavelet analysis

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @wave_clus_OpeningFcn, ...
    'gui_OutputFcn',  @wave_clus_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before wave_clus is made visible.
function wave_clus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for wave_clus
handles.output = hObject;
handles.datatype ='dg (joint)'; %dr change
handles.logf1=0;
handles.logf2=0;
set(handles.isi1_accept_button,'value',1);
set(handles.isi2_accept_button,'value',1);
set(handles.isi3_accept_button,'value',1);
set(handles.spike_shapes_button,'value',1);
set(handles.force_button,'value',0);
set(handles.plot_all_button,'value',1);
set(handles.plot_average_button,'value',0);
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
set(handles.feature1,'value',1);
set(handles.feature2,'value',2);
set(handles.feature3,'value',3);
set(handles.editSD,'value',1.5);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wave_clus wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wave_clus_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



clus_colors = [0 0 1; 1 0 0; 0 0.5 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
set(0,'DefaultAxesColorOrder',clus_colors)



% --- Executes on button press in load_data_button.
function load_data_button_Callback(hObject, eventdata, handles)
set(handles.isi1_accept_button,'value',1);
set(handles.isi2_accept_button,'value',1);
set(handles.isi3_accept_button,'value',1);
set(handles.isi1_reject_button,'value',0);
set(handles.isi2_reject_button,'value',0);
set(handles.isi3_reject_button,'value',0);
set(handles.isi1_nbins,'string','Auto');
set(handles.isi1_bin_step,'string','Auto');
set(handles.isi2_nbins,'string','Auto');
set(handles.isi2_bin_step,'string','Auto');
set(handles.isi3_nbins,'string','Auto');
set(handles.isi3_bin_step,'string','Auto');
set(handles.isi0_nbins,'string','Auto');
set(handles.isi0_bin_step,'string','Auto');
set(handles.force_button,'value',0);
set(handles.force_button,'string','Force');
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);

pack;
[filename, pathname] = uigetfile('joint*.mat','Select joint_ file');

if exist([pathname filename(7:end-3) 'dg'])==2
    handles.datatype='dg (joint)';
end

if exist([pathname filename(7:end-3) 'rhd'])==2
    handles.datatype='rhd (joint)';
end
 
switch char(handles.datatype)
    
    
    %dr change
    case 'dg (joint)'
       
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.org_pathname=pathname;
        handles.org_filename=filename;
        %handles.par = set_parameters_ascii_spikes(filename,handles);     %Load parameters
        handles.par.filename=filename;
        handles.par=set_parameters_master(handles.par);
        
        
        axes(handles.cont_data); cla
         
        %Load spikes and parameters
        eval(['load ' filename ';']);
        sub_directory=[pathname filename(7:end-4)];
        if exist('fls','var')==1
            file_to_cluster = char(fls(1));
        end
        load([file_to_cluster(1:end-3),'.mat']);
        
        handles.draq_p=params;
        handles.draq_d=data;
        handles.drta_p=drta_p;
        if exist('units_per_tet','var')
            handles.units_per_tet=units_per_tet;
        else
            for ii=1:4
                handles.units_per_tet(ii).no_units=0;
            end
        end
        handles.par.sr=params.srate;
        prompt = {'Enter channel or tetrode number:'};
        dlg_title = 'Input channel';
        num_lines = 1;
        def = {'1'};
        answr = inputdlg(prompt,dlg_title,num_lines,def);
        handles.drta_p.which_display=str2num(answr{1});
        which_display=str2num(answr{1});
        set(handles.tetrode,'string',answr{1});
        drta_p.which_display=handles.drta_p.which_display;
        spk_times=[];
        this_cluster_class=[];
        ii_spikes=1;
        if iscell(num_files)
            no_files=num_files{1};
        else
            no_files=num_files;
        end
        for filNum=1:no_files
            spk_times(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1)=all_timestamp(offset_chan_file(filNum,handles.drta_p.which_display)+1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)+noSpikesChFl(filNum,handles.drta_p.which_display));
            this_cluster_class(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1)=cluster_class(offset_chan_file(filNum,handles.drta_p.which_display) +1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)+noSpikesChFl(filNum,handles.drta_p.which_display));
            ii_spikes=ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display);
        end
%         spk_times=[];
%         spk_times(1:noSpikes(handles.drta_p.which_display))=all_timestamp(offset_for_chan(handles.drta_p.which_display)+...
%             1:offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
        spk_times=spk_times*1000;
        
        %Load clustering results
        cd(sub_directory)
        
        fname = [filename(7:end-4)];              %filename for interaction with SPC
        %clu=load(['times_ch' num2str(handles.drta_p.which_display) '_' fname '.dg_01.lab']);
        %Note: In this version clu is coming from times_ch
        load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.mat']);
         
        %The length of the dg_01 file was shortened because of a problem in
        %handling long file names by the c code in run_cluster
        if length(fname)>27
            tree=load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname(1:27) '.dg_01']);
        else
            tree=load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.dg_01']);
        end
        
        handles.par.fnamespc = fname;
        
        one_ch_file=['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.mat'];
        load(one_ch_file, 'spikes','inspk');
        pooled_spikes=[spikes(:,:,1) spikes(:,:,2) spikes(:,:,3) spikes(:,:,4)];
        handles.noSpikes=noSpikes;
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        USER_DATA{2} = pooled_spikes;
        nspk = size(pooled_spikes,1);
        if size(pooled_spikes,1)> handles.par.max_spk
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*handles.par.max_spk);
        else
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
        end
        set(handles.editSD,'string',num2str(handles.par.template_sdnum));
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        % Sets to zero fix buttons from aux figures
        for i=4:handles.par.max_clus
            eval(['handles.par.fix' num2str(i) '=0;'])
        end
        USER_DATA{1}=handles.par;
        USER_DATA{3} = spk_times(:)';
        USER_DATA{4} = clu;        USER_DATA{5} = tree;
        USER_DATA{6} = this_cluster_class(:)';
%         USER_DATA{6} = cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
%             offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA);
        cd (pathname);
        load([fname '.mat']);
        handles.draq_d=data;
        handles.draq_p=params;
        drta_p.which_display=handles.drta_p.which_display;
        handles.drta_p=drta_p;
        
        %Load jt_times to get the
        %This is done as try catch because with old Do_wave_clusdrdg
        %drta_p.f_strings was not saved
        try
            load(['jt_times_' file_to_cluster(1:end-3),'.mat']);
            set(handles.feature1,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
            set(handles.feature2,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
            set(handles.feature3,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
        catch
        end
        %Bring features chosen by Do_wave_clusdrdg
        %         for ii=1:10
        %            set(handles.feature1,'String')
        %         end
        
        
        % Update the handles structure
        guidata(hObject, handles);
        
    case 'rhd (joint)'

        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.org_pathname=pathname;
        handles.org_filename=filename;
        %handles.par = set_parameters_ascii_spikes(filename,handles);     %Load parameters
        handles.par.filename=filename;
        handles.par=set_parameters_master(handles.par);
        
        
        axes(handles.cont_data); cla
        
        %Load spikes and parameters
        eval(['load ' filename ';']);
        sub_directory=[pathname filename(7:end-4)];
        if exist('fls','var')==1
            file_to_cluster = char(fls(1));
        end
        load([file_to_cluster(1:end-4),'.mat']);
        handles.draq_p=params;
        handles.draq_d=data;
        handles.drta_p=drta_p;
        if exist('units_per_tet','var')
            handles.units_per_tet=units_per_tet;
        else
            for ii=1:4
                handles.units_per_tet(ii).no_units=0;
            end
        end
        handles.par.sr=params.srate;
        prompt = {'Enter channel or tetrode number:'};
        dlg_title = 'Input channel';
        num_lines = 1;
        def = {'1'};
        answr = inputdlg(prompt,dlg_title,num_lines,def);
        handles.drta_p.which_display=str2num(answr{1});
        which_display=str2num(answr{1});
        set(handles.tetrode,'string',answr{1});
        drta_p.which_display=handles.drta_p.which_display;
        spk_times=[];
        this_cluster_class=[];
        ii_spikes=1;
        for filNum=1:num_files
            spk_times(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1)=all_timestamp(offset_chan_file(filNum,handles.drta_p.which_display)+1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)+noSpikesChFl(filNum,handles.drta_p.which_display));
            this_cluster_class(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1)=cluster_class(offset_chan_file(filNum,handles.drta_p.which_display) +1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)+noSpikesChFl(filNum,handles.drta_p.which_display));
            ii_spikes=ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display);
        end
%         spk_times(1:noSpikes(handles.drta_p.which_display))=all_timestamp(offset_for_chan(handles.drta_p.which_display)+...
%             1:offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
         spk_times=spk_times*1000;
        
        %Load clustering results
        cd(sub_directory)
        
        fname = [filename(7:end-4)];              %filename for interaction with SPC
        %clu=load(['times_ch' num2str(handles.drta_p.which_display) '_' fname '.dg_01.lab']);
        %Note: In this version clu is coming from times_ch
        load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.mat']);
        
        %The length of the dg_01 file was shortened because of a problem in
        %handling long file names by the c code in run_cluster
        if length(fname)>28
            tree=load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname(1:28) '.dg_01']);
        else
            tree=load(['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.dg_01']);
        end
        
        handles.par.fnamespc = fname;
        
        one_ch_file=['times_tetr' num2str(handles.drta_p.which_display) '_' fname '.mat'];
        load(one_ch_file, 'spikes','inspk');
        pooled_spikes=[spikes(:,:,1) spikes(:,:,2) spikes(:,:,3) spikes(:,:,4)];
        handles.noSpikes=noSpikes;
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        USER_DATA{2} = pooled_spikes;
        nspk = size(pooled_spikes,1);
        if size(pooled_spikes,1)> handles.par.max_spk
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*handles.par.max_spk);
        else
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
        end
        set(handles.editSD,'string',num2str(handles.par.template_sdnum));
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        % Sets to zero fix buttons from aux figures
        for i=4:handles.par.max_clus
            eval(['handles.par.fix' num2str(i) '=0;'])
        end
        USER_DATA{1}=handles.par;
        USER_DATA{3} = spk_times(:)';
        USER_DATA{4} = clu;        USER_DATA{5} = tree;
%         USER_DATA{6} = cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
%             offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
        USER_DATA{6} = this_cluster_class(:)';
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA);
        cd (pathname);
        load([fname '.mat']);
        handles.draq_d=data;
        handles.draq_p=params;
        drta_p.which_display=handles.drta_p.which_display;
        handles.drta_p=drta_p;
        
        %Load jt_times to get the
        %This is done as try catch because with old Do_wave_clusdrdg
        %drta_p.f_strings was not saved
        try
            load(['jt_times_' file_to_cluster(1:end-3),'.mat']);
            set(handles.feature1,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
            set(handles.feature2,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
            set(handles.feature3,'String',{drta_p.f_strings{which_display,1},drta_p.f_strings{which_display,2},...
                drta_p.f_strings{which_display,3},drta_p.f_strings{which_display,4},drta_p.f_strings{which_display,5},...
                drta_p.f_strings{which_display,6},drta_p.f_strings{which_display,7},drta_p.f_strings{which_display,8},...
                drta_p.f_strings{which_display,9},drta_p.f_strings{which_display,10}});
        catch
        end
        %Bring features chosen by Do_wave_clusdrdg
        %         for ii=1:10
        %            set(handles.feature1,'String')
        %         end
        
        
        % Update the handles structure
        guidata(hObject, handles);
        
    case 'dra (joint)'
        [filename, pathname] = uigetfile('*.mat','Select joint_ file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        cd(pathname);
        handles.org_pathname=pathname;
        handles.org_filename=filename;
        %handles.par = set_parameters_ascii_spikes(filename,handles);     %Load parameters
        handles.par.filename=filename;
        handles.par=set_parameters_master(handles.par);
        
        axes(handles.cont_data); cla
        
        %Load spikes and parameters
        eval(['load ' filename ';']);
        sub_directory=[pathname filename(7:end-4)];
        if exist('fls','var')==1
            file_to_cluster = char(fls(1));
        end
        load([file_to_cluster,'.mat']);
        handles.draq_p=params;
        handles.draq_d=data;
        handles.drta_p=drta_p;
        handles.par.sr=params.srate;
        prompt = {'Enter channel number:'};
        dlg_title = 'Input channel';
        num_lines = 1;
        def = {'6'};
        answr = inputdlg(prompt,dlg_title,num_lines,def);
        handles.drta_p.which_display=str2num(answr{1});
        drta_p.which_display=handles.drta_p.which_display;
        spk_times=[];
        spk_times(1:noSpikes(handles.drta_p.which_display))=all_timestamp(offset_for_chan(handles.drta_p.which_display)+...
            1:offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
        spk_times=spk_times*1000;
        
        %Load clustering results
        cd(sub_directory)
        fname = [filename(7:end-4)];              %filename for interaction with SPC
        %clu=load(['times_ch' num2str(handles.drta_p.which_display) '_' fname '.dg_01.lab']);
        %Note: In this version clu is coming from times_ch
        load(['times_ch' num2str(handles.drta_p.which_display) '_' fname '.mat']);
        tree=load(['times_ch' num2str(handles.drta_p.which_display) '_' fname '.dg_01']);
        handles.par.fnamespc = fname;
        
        one_ch_file=['times_ch' num2str(handles.drta_p.which_display) '_' fname '.mat'];
        load(one_ch_file, 'spikes','inspk');
        handles.noSpikes=noSpikes;
        
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        USER_DATA{2} = spikes;
        nspk = size(spikes,1);
        if size(spikes,1)> handles.par.max_spk
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*handles.par.max_spk);
        else
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
        end
        set(handles.editSD,'string',num2str(handles.par.template_sdnum));
        set(handles.min_clus_edit,'string',num2str(handles.par.min_clus));
        % Sets to zero fix buttons from aux figures
        for i=4:handles.par.max_clus
            eval(['handles.par.fix' num2str(i) '=0;'])
        end
        
        USER_DATA{1}=handles.par;
        USER_DATA{3} = spk_times(:)';
        USER_DATA{4} = clu;        USER_DATA{5} = tree;
        USER_DATA{6} = cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
            offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display));
        USER_DATA{7} = inspk;
        set(handles.wave_clus_figure,'userdata',USER_DATA);
        cd (pathname);
        if isfield(handles.par,'correct_multichannel_artifact')
            if handles.par.correct_multichannel_artifact==1
                load([fname '.dra.mat']);
                %load([fname(1:end-4) '.dra.mat']);
            else
                load([fname '.dra.mat']);
            end
        else
            load([fname '.dra.mat']);
        end
        handles.draq_d=data;
        handles.draq_p=params;
        drta_p.which_display=handles.drta_p.which_display;
        handles.drta_p=drta_p;
        
        
        % Update the handles structure
        guidata(hObject, handles);
        
end


temp=find_temp(tree,handles);                                   %Selects temperature.
temperature=handles.par.mintemp+temp*handles.par.tempstep;
axes(handles.temperature_plot);
switch handles.par.temp_plot
    case 'lin'
        plot([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [handles.par.min_clus handles.par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
    case 'log'
        semilogy([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [handles.par.min_clus handles.par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
end
xlim([0 handles.par.maxtemp])
xlabel('Temperature');
if handles.par.temp_plot == 'log'
    set(get(gca,'ylabel'),'vertical','Cap');
else
    set(get(gca,'ylabel'),'vertical','Baseline');
end
ylabel('Clusters size');
set(handles.file_name,'string',[pathname filename]);

if size(clu,2)-2 < size(spikes,1);
    classes = clu(temp,3:end)+1;
    classes = [classes(:)' zeros(1,size(spikes,1)-handles.par.max_spk)];
else
    classes = clu(temp,3:end)+1;
end

guidata(hObject, handles);
USER_DATA = get(handles.wave_clus_figure,'userdata');
USER_DATA{6} = classes(:)';
USER_DATA{8} = temp;
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.
set(handles.wave_clus_figure,'userdata',USER_DATA);

handles.setclus = 0;
plot_spikes(handles);


% --- Executes on button press in change_temperature_button.
function change_temperature_button_Callback(hObject, eventdata, handles)
axes(handles.temperature_plot)
[temp, aux]= ginput(1);                                          %gets the mouse input
temp = round((temp-handles.par.mintemp)/handles.par.tempstep);
if temp < 1; temp=1;end                                         %temp should be within the limits
if temp > handles.par.num_temp; temp=handles.par.num_temp; end
min_clus = round(aux);
set(handles.min_clus_edit,'string',num2str(min_clus));

USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.min_clus = min_clus;
clu = USER_DATA{4};
classes = clu(temp,3:end)+1;
tree = USER_DATA{5};
USER_DATA{1} = par;
USER_DATA{6} = classes(:)';
USER_DATA{8} = temp;
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.
set(handles.wave_clus_figure,'userdata',USER_DATA);
temperature=handles.par.mintemp+temp*handles.par.tempstep;

switch par.temp_plot
    case 'lin'
        plot([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep],[par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
    case 'log'
        semilogy([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
end
xlim([0 handles.par.maxtemp])
xlabel('Temperature');
if par.temp_plot == 'log'
    set(get(gca,'ylabel'),'vertical','Cap');
else
    set(get(gca,'ylabel'),'vertical','Baseline');
end
ylabel('Clusters size');
handles.setclus = 0;
plot_spikes(handles);
set(handles.force_button,'value',0);
set(handles.force_button,'string','Force');
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end


% --- Change min_clus_edit
function min_clus_edit_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.min_clus = str2num(get(hObject, 'String'));
clu = USER_DATA{4};
temp = USER_DATA{8};
classes = clu(temp,3:end)+1;
tree = USER_DATA{5};
USER_DATA{1} = par;
USER_DATA{6} = classes(:)';
USER_DATA{9} = classes(:)';                                     %backup for non-forced classes.
set(handles.wave_clus_figure,'userdata',USER_DATA);
temperature=handles.par.mintemp+temp*handles.par.tempstep;

axes(handles.temperature_plot)
switch par.temp_plot
    case 'lin'
        plot([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep],[par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
    case 'log'
        semilogy([handles.par.mintemp handles.par.maxtemp-handles.par.tempstep], ...
            [par.min_clus par.min_clus],'k:',...
            handles.par.mintemp+(1:handles.par.num_temp)*handles.par.tempstep, ...
            tree(1:handles.par.num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
end
xlim([0 handles.par.maxtemp])
xlabel('Temperature');
if par.temp_plot == 'log'
    set(get(gca,'ylabel'),'vertical','Cap');
else
    set(get(gca,'ylabel'),'vertical','Baseline');
end
ylabel('Clusters size');
handles.setclus = 0;
plot_spikes(handles);
set(handles.force_button,'value',0);
set(handles.force_button,'string','Force');
set(handles.fix1_button,'value',0);
set(handles.fix2_button,'value',0);
set(handles.fix3_button,'value',0);
for i=4:par.max_clus
    eval(['par.fix' num2str(i) '=0;']);
end



% --- Executes on button press in save_clusters_button.
function save_clusters_button_Callback(hObject, eventdata, handles)

cd(handles.org_pathname);
set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
set(gcf,'paperposition',[.25 .25 10.5 7.8])
eval(['print -djpeg40 waveclus_tetr_' num2str(handles.drta_p.which_display) '_' handles.org_filename '.jpg']);

switch char(handles.datatype)
    case 'dg (joint)'
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        classes = USER_DATA{6};
        
        cont=0;
        
        % Classes should be consecutive numbers
        i=1;
        while i<=max(classes)
            if isempty(classes(find(classes==i)))
                for k=i+1:max(classes)
                    classes(find(classes==k))=k-1;
                end
            else
                i=i+1;
            end
        end
        
        
        %Retrieve all the original data from file
        cd(handles.org_pathname);
        load(handles.org_filename);
        chNo=handles.drta_p.which_display;
        
        %         cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
        %             offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display))=classes';
        %
        ii_spikes=1;
        if iscell(num_files)
            no_files=num_files{1};
        else
            no_files=num_files;
        end
        for filNum=1:no_files
            cluster_class(offset_chan_file(filNum,handles.drta_p.which_display)+1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)...
                +noSpikesChFl(filNum,handles.drta_p.which_display))=...
                classes(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1);
            ii_spikes=ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display);
        end
        
        units_per_tet=handles.units_per_tet;
        
        %save joint file
        if exist('fls','var')==1
            save(handles.org_filename, 'num_files', 'fls','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes','units_per_tet');
            no_files=num_files{1};
        else
            save(handles.org_filename, 'num_files', 'file_to_cluster','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes','units_per_tet');
             no_files=num_files;
        end
        
        noSpikes=[];
        for filNum=1:no_files
            
            if exist('fls','var')==1
                dg_file=char(fls(filNum));
                jt_times_file=['jt_times_' dg_file(1:end-3) '.mat'];
            else
                jt_times_file=['jt_times_' file_to_cluster(1:end-3) '.mat'];
            end
            
            
            load(jt_times_file);
            noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
            all_timestamp_per_file=[];
            cluster_class_per_file=[];
            
            for jj=1:4
                szat=size(all_timestamp_per_file);
                offset_for_chan(jj)=szat(2);
                all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
                cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
            end
            
            
            if isfield(drta_p,'tetr_processed')==0
                drta_p.tetr_processed=[0 0 0 0];
            end
            drta_p.tetr_processed(chNo)=1;
            
            units_per_tet=handles.units_per_tet;
            
            %The user may have changed the directroy name (to fix the
            %annoying space problem)
            drta_p.PathName=[directory '/'];
            drta_p.fullName=[drta_p.PathName drta_p.FileName];
            
            if isfield(handles.par,'doBehavior')
                
                if handles.par.doBehavior==1
                    %                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','lickbit','dtime');
                    %                else
                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','units_per_tet');
                end
                
            else
                save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','units_per_tet');
            end
            
            %Now generate the drg file
            drgRead_jt_times(handles.org_pathname,jt_times_file)
        end
        handles.org_filename
    
    case 'rhd (joint)'
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        classes = USER_DATA{6};
        
        cont=0;
        
        % Classes should be consecutive numbers
        i=1;
        while i<=max(classes)
            if isempty(classes(find(classes==i)))
                for k=i+1:max(classes)
                    classes(find(classes==k))=k-1;
                end
            else
                i=i+1;
            end
        end
        
        
        %Retrieve all the original data from file
        cd(handles.org_pathname);
        load(handles.org_filename);
        chNo=handles.drta_p.which_display;
        
        ii_spikes=1;
        for filNum=1:num_files
            cluster_class(offset_chan_file(filNum,handles.drta_p.which_display)+1 ...
                :offset_chan_file(filNum,handles.drta_p.which_display)...
                +noSpikesChFl(filNum,handles.drta_p.which_display))=...
                classes(ii_spikes:ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display)-1);
            ii_spikes=ii_spikes+noSpikesChFl(filNum,handles.drta_p.which_display);
        end
        
%         cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
%             offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display))=classes';
%         
        units_per_tet=handles.units_per_tet;
        
        %save joint file
        if exist('fls','var')==1
            save(handles.org_filename, 'fls','num_files','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes','units_per_tet');
        else
            save(handles.org_filename, 'num_files', 'file_to_cluster','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes','units_per_tet');
        end
        
        noSpikes=[];
        for filNum=1:num_files
            
            if exist('fls','var')==1
                dg_file=fls{filNum};
                jt_times_file=['jt_times_' dg_file(1:end-4) '.mat'];
            else
                jt_times_file=['jt_times_' file_to_cluster(1:end-4) '.mat'];
            end
             
            
            load(jt_times_file);
            noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
            all_timestamp_per_file=[];
            cluster_class_per_file=[];
            
            for jj=1:4
                szat=size(all_timestamp_per_file);
                offset_for_chan(jj)=szat(2);
                all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
                cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
            end
            
            
            if isfield(drta_p,'tetr_processed')==0
                drta_p.tetr_processed=[0 0 0 0];
            end
            drta_p.tetr_processed(chNo)=1;
            
            units_per_tet=handles.units_per_tet;
             
            %The user may have changed the directroy name (to fix the
            %annoying space problem)
            drta_p.PathName=[directory '/'];
            drta_p.fullName=[drta_p.PathName drta_p.FileName];
            
            if isfield(handles.par,'doBehavior')
                
                if handles.par.doBehavior==1
                    %                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','lickbit','dtime');
                    %                else
                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','units_per_tet');
                end
                
            else
                save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','units_per_tet');
            end
            
            %Now generate the drg file
            drgRead_jt_times(handles.org_pathname,jt_times_file)
        end
        handles.org_filename
        
    case 'dra (joint)'
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        
        classes = USER_DATA{6};
        
        cont=0;
        
        % Classes should be consecutive numbers
        i=1;
        while i<=max(classes)
            if isempty(classes(find(classes==i)))
                for k=i+1:max(classes)
                    classes(find(classes==k))=k-1;
                end
            else
                i=i+1;
            end
        end
        
        
        %Retrieve all the original data from file
        cd(handles.org_pathname);
        load(handles.org_filename);
        chNo=handles.drta_p.which_display;
        
        cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
            offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display))=classes';
        if exist('fls','var')==1
            save(handles.org_filename, 'num_files', 'fls','directory','sub_directory','offset_chan_file', 'noSpikesChFl','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes');
        else
            save(handles.org_filename, 'num_files', 'file_to_cluster','directory','sub_directory','offset_chan_file', 'noSpikesChFl','all_timestamp',...
                'cluster_class','offset_for_chan', 'noSpikes');
        end
        
        noSpikes=[];
        for filNum=1:num_files{1}
            
            if exist('fls','var')==1
                dra_file=char(fls(filNum));
                
                jt_times_file=['jt_times_' dra_file(1:end-4) '.mat'];
            else
                jt_times_file=['jt_times_' file_to_cluster(1:end-4) '.mat'];
            end
            
            
            load(jt_times_file);
            noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
            all_timestamp_per_file=[];
            cluster_class_per_file=[];
            
            for jj=1:4
                szat=size(all_timestamp_per_file);
                offset_for_chan(jj)=szat(2);
                all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
                cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
            end
            
            if isfield(handles.par,'doBehavior')
                
                if handles.par.doBehavior==1
                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','lickbit','dtime');
                else
                    save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d');
                end
                
            else
                save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d');
            end
            
        end
        handles.org_filename
        
end

msgbox('Saved jt_times for this tetrode')



%SETTING OF FORCE MEMBERSHIP
% --------------------------------------------------------------------
function force_button_Callback(hObject, eventdata, handles)
%set(gcbo,'value',1);
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
spikes = USER_DATA{2};
classes = USER_DATA{6};
inspk = USER_DATA{7};

% Fixed clusters are not considered for forcing
if get(handles.fix1_button,'value') ==1
    fix_class = USER_DATA{10}';
    classes(fix_class)=-1;
end
if get(handles.fix2_button,'value') ==1
    fix_class = USER_DATA{11}';
    classes(fix_class)=-1;
end
if get(handles.fix3_button,'value') ==1
    fix_class = USER_DATA{12}';
    classes(fix_class)=-1;
end
% Get fixed clusters from aux figures
for i=4:par.max_clus
    eval(['fixx = par.fix' num2str(i) ';']);
    if fixx == 1
        fix_class = USER_DATA{12+i-3}';
        classes(fix_class)=-i;
    end
end


switch handles.par.force_feature
    case 'spk'
        f_in  = spikes(find(classes~=0 & classes~=-1),:);
        f_out = spikes(find(classes~=-1),:);
    case 'wav'
        if isempty(inspk)
            [inspk, handles.drta_p.f_strings] = wave_features_dr(spikes,handles);        % Extract spike features.
            USER_DATA{7} = inspk;
        end
        f_in  = inspk(find(classes~=0 & classes~=-1),:);
        f_out = inspk(find(classes~=-1),:);
end
class_in = classes(find(classes~=0 & classes~=-1));

if get(handles.force_button,'value') ==1
    class_out = force_membership_wc(f_in, class_in, f_out, handles);
    classes(find(classes~=-1)) = class_out;
    set(handles.force_button,'string','Forced');
elseif get(handles.force_button,'value') ==0
    classes = USER_DATA{9};
    set(handles.force_button,'string','Force');
end
USER_DATA{6} = classes(:)';
set(handles.wave_clus_figure,'userdata',USER_DATA)

handles.setclus = 1;
plot_spikes(handles);

pfft=1

% set(handles.fix1_button,'value',0);
% set(handles.fix2_button,'value',0);
% set(handles.fix3_button,'value',0);
% for i=4:par.max_clus
%     eval(['par.fix' num2str(i) '=0;']);
% end




% PLOT ALL PROJECTIONS BUTTON
% --------------------------------------------------------------------
function Plot_all_projections_button_Callback(hObject, eventdata, handles)
Plot_all_features(handles)
% --------------------------------------------------------------------


% fix1 button --------------------------------------------------------------------
function fix1_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==1);
if get(handles.fix1_button,'value') ==1
    USER_DATA{10} = fix_class;
else
    USER_DATA{10} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


% fix2 button --------------------------------------------------------------------
function fix2_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==2);
if get(handles.fix2_button,'value') ==1
    USER_DATA{11} = fix_class;
else
    USER_DATA{11} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


% fix3 button --------------------------------------------------------------------
function fix3_button_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
fix_class = find(classes==3);
if get(handles.fix3_button,'value') ==1
    USER_DATA{12} = fix_class;
else
    USER_DATA{12} = [];
end
set(handles.wave_clus_figure,'userdata',USER_DATA)
h_figs=get(0,'children');
h_fig2 = findobj(h_figs,'tag','wave_clus_aux');
h_fig1 = findobj(h_figs,'tag','wave_clus_aux1');
set(h_fig2,'userdata',USER_DATA)
set(h_fig1,'userdata',USER_DATA)


%SETTING OF SPIKE FEATURES OR PROJECTIONS
% --------------------------------------------------------------------
function spike_shapes_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.spike_features_button,'value',0);
handles.setclus = 1;
plot_spikes(handles);
% -------------------------------------------------------------------
function spike_features_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.spike_shapes_button,'value',0);
handles.setclus = 1;
plot_spikes(handles);


%SETTING OF SPIKE PLOTS
% --------------------------------------------------------------------
function plot_all_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.plot_average_button,'value',0);
handles.setclus = 1;
plot_spikes(handles);
% --------------------------------------------------------------------
function plot_average_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.plot_all_button,'value',0);
handles.setclus = 1;
plot_spikes(handles);



%SETTING OF ISI HISTOGRAMS
% --------------------------------------------------------------------
function isi1_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.nbins1 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi1_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.bin_step1 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi2_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.nbins2 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi2_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.bin_step2 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi3_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.nbins3 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi3_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.bin_step3 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi0_nbins_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.nbins0 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)
% --------------------------------------------------------------------
function isi0_bin_step_Callback(hObject, eventdata, handles)
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
par.bin_step0 = str2num(get(hObject, 'String'));
USER_DATA{1} = par;
set(handles.wave_clus_figure,'userdata',USER_DATA);
handles.setclus = 1;
plot_spikes(handles)



%SETTING OF ISI BUTTONS

% --------------------------------------------------------------------
function isi1_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi1_reject_button,'value',0);

% --------------------------------------------------------------------
function isi1_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi1_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
classes(find(classes==1))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
set(handles.wave_clus_figure,'userdata',USER_DATA);

handles.setclus = 1;
plot_spikes(handles)

set(gcbo,'value',0);
set(handles.isi1_accept_button,'value',1);

% --------------------------------------------------------------------
function isi2_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi2_reject_button,'value',0);

% --------------------------------------------------------------------
function isi2_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi2_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
classes(find(classes==2))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
set(handles.wave_clus_figure,'userdata',USER_DATA);

handles.setclus = 1;
plot_spikes(handles)

set(gcbo,'value',0);
set(handles.isi2_accept_button,'value',1);

% --------------------------------------------------------------------
function isi3_accept_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi3_reject_button,'value',0);

% --------------------------------------------------------------------
function isi3_reject_button_Callback(hObject, eventdata, handles)
set(gcbo,'value',1);
set(handles.isi3_accept_button,'value',0);
USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};
classes(find(classes==3))=0;
USER_DATA{6} = classes;
USER_DATA{9} = classes;
set(handles.wave_clus_figure,'userdata',USER_DATA);

handles.setclus = 1;
plot_spikes(handles)

set(gcbo,'value',0);
set(handles.isi3_accept_button,'value',1);


% --- Executes during object creation, after setting all properties.
function isi1_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi1_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi1_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi1_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi2_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi2_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi2_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi2_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi3_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi3_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi3_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi3_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi0_nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi0_nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function isi0_bin_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi0_bin_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function min_clus_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_clus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature1.
function feature1_Callback(hObject, eventdata, handles)
% hObject    handle to feature1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns feature1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature1
handles.setclus = 1;
plot_spikes(handles);

% --- Executes during object creation, after setting all properties.
function feature1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature2.
function feature2_Callback(hObject, eventdata, handles)
% hObject    handle to feature2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns feature2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature2
handles.setclus = 1;
plot_spikes(handles);

% --- Executes during object creation, after setting all properties.
function feature2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature3.
function feature3_Callback(hObject, eventdata, handles)
% hObject    handle to feature3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns feature3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature3


% --- Executes during object creation, after setting all properties.
function feature3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feature3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in threeD.
function threeD_Callback(hObject, eventdata, handles)
% hObject    handle to threeD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(11);
USER_DATA = get(handles.wave_clus_figure,'userdata');
inspk = USER_DATA{7};
clusters=USER_DATA{6};
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];
hold on
for clu=min(clusters):max(clusters)
    this_clu=find(clusters==clu);
    plot3(inspk(this_clu,get(handles.feature1,'value')),inspk(this_clu,get(handles.feature2,'value')),inspk(this_clu,get(handles.feature3,'value')),['.' colors(clu+1)],'MarkerSize',2);
    xstr=get(handles.feature1,'string');
    xlabel(xstr(get(handles.feature1,'value')));
    ystr=get(handles.feature2,'string');
    ylabel(ystr(get(handles.feature2,'value')));
    zstr=get(handles.feature3,'string');
    zlabel(zstr(get(handles.feature3,'value')));
end
hold off



function editSD_Callback(hObject, eventdata, handles)
% hObject    handle to editSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSD as text
%        str2double(get(hObject,'String')) returns contents of editSD as a double
handles.par.template_sdnum=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in join_clusters.
function join_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to join_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Cluster:','Cluster:','Cluster:','Cluster:'};
dlg_title = 'Enter clusters to be joined (-1=no entry)';
num_lines = 1;
def = {'-1','-1','-1','-1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

USER_DATA = get(handles.wave_clus_figure,'userdata');
classes = USER_DATA{6};

for ii=1:4
    num_answer(ii)=str2num(answer{ii});
end
clus_tojoin=find(num_answer~=-1);
new_clus=min(num_answer(clus_tojoin));
for ii=1:4
    is_there=find(clus_tojoin==ii)
    if ~isempty(is_there)
        classes(find(classes==num_answer(ii)))=new_clus;
    end
end

USER_DATA{6} = classes(:);
set(handles.wave_clus_figure,'userdata',USER_DATA)

handles.setclus = 0;
plot_spikes(handles);


% --- Executes on button press in stemPlot.
function stemPlot_Callback(hObject, eventdata, handles)
% hObject    handle to stemPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    close(6)
catch
end
figure(6)


colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];



%Retrieve all the original data from file
cd(handles.org_pathname);
load(handles.org_filename);
chNo=handles.drta_p.which_display;

USER_DATA = get(handles.wave_clus_figure,'userdata');

classes = USER_DATA{6};

cluster_class(offset_for_chan(handles.drta_p.which_display)+1:...
    offset_for_chan(handles.drta_p.which_display)+noSpikes(handles.drta_p.which_display))=classes';


noSpikes=[];
ylims=[];
for filNum=1:num_files{1}(1)
    
    
    noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
    all_timestamp_per_file=[];
    cluster_class_per_file=[];
    
    
    szat=size(all_timestamp_per_file);
    offset_for_chan(chNo)=szat(2);
    all_timestamp_per_file=[];
    all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,chNo)+1:offset_chan_file(filNum,chNo)+noSpikes(chNo))];
    cluster_class_per_file=[];
    cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,chNo)+1:offset_chan_file(filNum,chNo)+noSpikes(chNo))];
    
    classNo=max(cluster_class_per_file)
    
    
    for cNo=1: classNo
        subplot(double(classNo),double(num_files{1}(1)),double(filNum+double(num_files{1}(1))*(cNo-1)))
        try
        histogram(all_timestamp_per_file(cluster_class_per_file==cNo),...
            [min(all_timestamp_per_file(cluster_class_per_file==cNo)):0.1:max(all_timestamp_per_file(cluster_class_per_file==cNo))],...
            'FaceColor',colors(cNo+1),'EdgeColor',colors(cNo+1));
        yl=ylim;
        ylims(cNo,filNum)=yl(2);
        catch
        end
        title(['file' num2str(filNum)])
    end
    
end

for filNum=1:num_files{1}(1)
    for cNo=1: classNo
        subplot(double(classNo),double(num_files{1}(1)),double(filNum+double(num_files{1}(1))*(cNo-1)))
        ylim([0 max(ylims(cNo,:))])
    end
end
    


% --- Executes on button press in enterROI.
function enterROI_Callback(hObject, eventdata, handles)
% hObject    handle to enterROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=msgbox('Now re-size plot by entering two points, then enter ROI points with a return at end','ROI');

try
    close 9
catch
end

h=figure(9);
USER_DATA = get(handles.wave_clus_figure,'userdata');
inspk = USER_DATA{7};
clusters=USER_DATA{6};
last_cluster=max(clusters);
colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];
hold on
for clu=min(clusters):max(clusters)
    this_clu=find(clusters==clu);
    %     plot(log10(inspk(this_clu,get(handles.feature1,'value'))),...
    %         log10(inspk(this_clu,get(handles.feature2,'value'))),['.' colors(clu+1)],'MarkerSize',2);
    %     plot(inspk(this_clu,get(handles.feature1,'value')),...
    %         inspk(this_clu,get(handles.feature2,'value')),['.' colors(clu+1)],'MarkerSize',2);
    if (handles.logf1==0)&(handles.logf2==0)
        plot(inspk(this_clu,get(handles.feature1,'value')),...
            inspk(this_clu,get(handles.feature2,'value')),['.' colors(clu+1)],'MarkerSize',1);
    end
    if (handles.logf1==1)&(handles.logf2==0)
        plot(log10(inspk(this_clu,get(handles.feature1,'value'))),...
            inspk(this_clu,get(handles.feature2,'value')),['.' colors(clu+1)],'MarkerSize',1);
    end
    if (handles.logf1==0)&(handles.logf2==1)
        plot(inspk(this_clu,get(handles.feature1,'value')),...
            log10(inspk(this_clu,get(handles.feature2,'value'))),['.' colors(clu+1)],'MarkerSize',1);
    end
    if (handles.logf1==1)&(handles.logf2==1)
        plot(log10(inspk(this_clu,get(handles.feature1,'value'))),...
            log10(inspk(this_clu,get(handles.feature2,'value'))),['.' colors(clu+1)],'MarkerSize',1);
    end
    xstr=get(handles.feature1,'string');
    xlabel(xstr(get(handles.feature1,'value')));
    ystr=get(handles.feature2,'string');
    ylabel(ystr(get(handles.feature2,'value')));
end


%First re-size the plot
for ii=1:2
    [f1(ii), f2(ii)]= ginput(1);
end

if f1(1)>f1(2)
    xlim([f1(2) f1(1)])
else
    xlim([f1(1) f1(2)])
end

if f2(1)>f2(2)
    ylim([f2(2) f2(1)])
else
    ylim([f2(1) f2(2)])
end

%Then get the ROI
ii=1;
[f1, f2]= ginput(1);
f1_vec(ii)=f1;
f2_vec(ii)=f2;
while ~isempty(f1)
    [f1, f2]= ginput(1);
    if ~isempty(f1)
        ii=ii+1;
        f1_vec(ii)=f1;
        f2_vec(ii)=f2;
        line(f1_vec(ii-1:ii),f2_vec(ii-1:ii))
    end
    
end
line([f1_vec(ii) f1_vec(1)],[f2_vec(ii),f2_vec(1)])

%I had a problem with inploygon missing a few points
f1_vec(end+1)=f1_vec(1);
f1_vec=[f1_vec f1_vec];
f2_vec(end+1)=f2_vec(1);
f2_vec=[f2_vec f2_vec];

f1length=length(f1_vec);
maxf1=max(f1_vec);
minf1=min(f1_vec);
maxf2=max(f2_vec);
minf2=min(f2_vec);

%Get user_data
USER_DATA = get(handles.wave_clus_figure,'userdata');
inspk = USER_DATA{7};
szinspk=size(inspk);

%Set clusters
if (handles.logf1==0)&(handles.logf2==0)
    clusters(inpolygon(inspk(:,get(handles.feature1,'value')),inspk(:,get(handles.feature2,'value')),f1_vec,f2_vec))=last_cluster+1;
end
if (handles.logf1==1)&(handles.logf2==0)
    clusters(inpolygon(log10(inspk(:,get(handles.feature1,'value'))),inspk(:,get(handles.feature2,'value')),f1_vec,f2_vec))=last_cluster+1;
end
if (handles.logf1==0)&(handles.logf2==1)
    clusters(inpolygon(inspk(:,get(handles.feature1,'value')),log10(inspk(:,get(handles.feature2,'value'))),f1_vec,f2_vec))=last_cluster+1;
end
if (handles.logf1==1)&(handles.logf2==1)
    clusters(inpolygon(log10(inspk(:,get(handles.feature1,'value'))),log10(inspk(:,get(handles.feature2,'value'))),f1_vec,f2_vec))=last_cluster+1;
end

USER_DATA{6}=clusters;
USER_DATA{9} = clusters;


set(handles.wave_clus_figure,'userdata',USER_DATA)

handles.setclus = 1;
plot_spikes(handles);


% --- Executes on button press in logf1.
function logf1_Callback(hObject, eventdata, handles)
% hObject    handle to logf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logf1

handles.logf1=get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in logf2.
function logf2_Callback(hObject, eventdata, handles)
% hObject    handle to logf2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logf2

handles.logf2=get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in resize.
function resize_Callback(hObject, eventdata, handles)
% hObject    handle to resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resize


% --- Executes on selection change in force_feature.
function force_feature_Callback(hObject, eventdata, handles)
% hObject    handle to force_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns force_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from force_feature

if get(hObject,'Value')==1
    handles.par.force_feature='spk';
else
    handles.par.force_feature='wav';
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function force_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
