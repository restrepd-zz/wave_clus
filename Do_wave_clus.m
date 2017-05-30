function varargout = Do_wave_clus(varargin)
% DO_WAVE_CLUS MATLAB code for Do_wave_clus.fig
%      DO_WAVE_CLUS, by itself, creates a new DO_WAVE_CLUS or raises the existing
%      singleton*.
%
%      H = DO_WAVE_CLUS returns the handle to a new DO_WAVE_CLUS or the handle to
%      the existing singleton*.
%
%      DO_WAVE_CLUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DO_WAVE_CLUS.M with the given input arguments.
%
%      DO_WAVE_CLUS('Property','Value',...) creates a new DO_WAVE_CLUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Do_wave_clus_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Do_wave_clus_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% This is a modification of the wave_clus code from Quian Quiroga
% for use in the Restrepo lab
%
% The Quian Quiroga code is available in
%
% http://www2.le.ac.uk/centres/csn/research-2/spike-sorting
%
% We have added the ablity to use peak to valley and PCA in addition to the
% wavelet analysis
%

% Edit the above text to modify the response to help Do_wave_clus

% Last Modified by GUIDE v2.5 28-May-2017 17:05:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Do_wave_clus_OpeningFcn, ...
    'gui_OutputFcn',  @Do_wave_clus_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Do_wave_clus is made visible.
function Do_wave_clus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Do_wave_clus (see VARARGIN)

% Choose default command line output for Do_wave_clus
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Do_wave_clus wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Do_wave_clus_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_run_dg.
function load_run_dg_Callback(hObject, eventdata, handles)
% hObject    handle to load_run_dg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% %GET SPIKES prameters
tic

try
    close 1
catch
end

try
    close 2
catch
end
handles.file_name='';
handles.datatype ='CSC data'; %Needed to avoid output to figure
handles.flag=0;

print2file =1;                              %for saving printouts.
handles.par.prog = 'do_wave_clus';
handles.par=set_parameters_master(handles.par);
handles.par.dgorrhd=1;
directory=handles.directory;

bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');
pos1  = [bdwidth+(1/6)*scnsize(4),...
    (3/8)*scnsize(4),...
    scnsize(3)*(2/3) - 2*bdwidth,...
    scnsize(4)*(1/2) - (topbdwidth + bdwidth)];


set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8])

cd(directory);
num_files=handles.numfiles;



cd(directory);


handles.drta_p.draName=handles.files{1};
handles.drta_p.fullName=handles.files{1};
file_to_cluster = handles.files{1};
%     end
if ismac
    sub_directory=[directory '/' handles.drta_p.fullName(1:end-3)];
else
    sub_directory=[directory '\' handles.drta_p.fullName(1:end-3)];
end
load([handles.drta_p.draName(1:end-3),'.mat']);

an_error=0;
try
    handles.draq_p=params;
catch
    an_error=1;
end

if an_error==0
    handles.drta_p=drta_p;
    
    noSpikes=[];
    spikesExcluded=[];
    timestamp_offset=0;
    all_timestamp=[];   %Timestamps
    offset_for_chan=[]; %For each channel this has the offset within all_timestamp
    cluster_class=[];
    
    offset_chan_file=zeros(handles.numfiles,4);
    
    noSpikesChFl=zeros(handles.numfiles,4);
    noSpikesChFlExc=zeros(handles.numfiles,4);
    f_strings=cell(4,10);
    
    %Process data for each tetrode
    for tets=1:4
        tetrode_number=tets
        %For each tetrode find spikes, spike features and do processing
        
        spikes=[];
        index=[];
        timestamp=[];
        nexc=0;
        lbl='detect'
        
        
        
        for filNum=1:handles.numfiles
            file_number=filNum
            noSpikesChFlExc(filNum,tets)=0;
            handles.drta_p.draName=handles.files{filNum};
            cd(directory);
            load([handles.drta_p.draName(1:end-3),'.mat']);
            handles.draq_p=params;
            handles.par.sr=handles.draq_p.ActualRate;
            if handles.par.read_entire_file==1
                handles.draq_d=data;
                handles.draq_d.data=wc_read_dg(handles);
            end
            handles.drta_p=drta_p;
            handles.drta_p.draName=handles.files{filNum};
            handles.drta_p.fullName=handles.files{filNum};
            handles.par.filename=handles.files{filNum};
            handles.drta_p.tets=tets;
            sz_ts=size(timestamp);
            file_start_sz=sz_ts(2);
            
            for ii=1:handles.draq_d.noTrials
                %for ii=1:2
                trial_number=ii
                
                
                
                handles.drta_p.trialNo=ii;
                all_tet_ch_proc=1;
                for chno=4*(tets-1)+1:4*(tets-1)+4
                    
                    if (handles.drta_p.trial_ch_processed(chno,handles.drta_p.trialNo)==0)
                        all_tet_ch_proc=0;
                    end
                    
                end
                
                
                if (all_tet_ch_proc==1)
                    %Do not process if one of the channels for the tetrode was
                    %taken off for this trial
                    
                    [aux_spikes,thr,aux_index,aux_excluded,aux_nspk] = drdg_amp_detect_wc(handles);     % Detection with amp. thresh.
                    
                    if aux_nspk>0
                        shift=handles.draq_d.t_trial(ii)*handles.draq_p.ActualRate;
                        spikes=[spikes;aux_spikes];
                        index=[index aux_index];
                        timestamp=[timestamp (shift+aux_index)/handles.draq_p.ActualRate];
                        nexc=nexc+aux_excluded;
                        noSpikesChFlExc(filNum,tets)=noSpikesChFlExc(filNum,tets)+aux_excluded;
                    end
                end
                
                
            end
            sz_ts=size(timestamp);
            file_end_sz=sz_ts(2);
            offset_chan_file(filNum,tets)=timestamp_offset+file_start_sz;
            noSpikesChFl(filNum,tets)=file_end_sz-file_start_sz;
            %fclose(handles.drta_p.fullName);
        end
        
        
        
        nspk = size(spikes,1);
        %feature('memstats')
        all_timestamp(timestamp_offset+1:timestamp_offset+nspk)=timestamp(:);
        offset_for_chan(tets)=timestamp_offset;
        noSpikes(tets)=nspk;
        spikesExcluded(tets)=nexc;
        
        %Cluster starts here
        mkdir(sub_directory);
        cd(sub_directory);
        lbl='cluster'
        
        timestamp=timestamp*1000;
        handles.par.fname = ['data_' file_to_cluster];   %filename for interaction with SPC
        
        if nspk > handles.par.max_spk
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*handles.par.max_spk);
        else
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
        end
        
        % CALCULATES INPUTS TO THE CLUSTERING ALGORITHM.
        if nspk>0
            [inspk, f_str] = wave_features_drdg(spikes,handles);              %takes wavelet coefficients.
            for ii=1:10
                f_strings{tets,ii}=f_str{ii};
            end
            
            % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
            if size(spikes,1)> handles.par.max_spk
                naux = floor(handles.par.max_spk/handles.numfiles);
                inspk_aux=[];
                for filNum=1:handles.numfiles
                    if naux>noSpikesChFl(filNum,tets)
                        naux2(filNum)=noSpikesChFl(filNum,tets);
                    else
                        naux2(filNum)=naux;
                    end
                    szin=size(inspk_aux);
                    offset_aux(filNum)=szin(1);
                    inspk_aux = [inspk_aux; inspk(1+offset_chan_file(filNum,tets)-offset_chan_file(1,tets):...
                        offset_chan_file(filNum,tets)-offset_chan_file(1,tets)+naux2(filNum),:)];
                    
                end
            else
                inspk_aux = inspk;
            end
            
            
            %INTERACTION WITH SPC
            save(handles.par.fname_in,'inspk_aux','-ascii');
            [aux_clu, tree] = run_cluster(handles);
            [temp] = find_temp(tree,handles);
            
            
            %DEFINE CLUSTERS
            szclu=size(aux_clu);
            szspikes=size(spikes);
            clu=-1*ones(szclu(1),szspikes(1)+2);
            clu(:,1:2)=aux_clu(:,1:2);
            
            if size(spikes,1)> handles.par.max_spk
                for filNum=1:handles.numfiles
                    
                    clu(:,offset_chan_file(filNum,tets)+1-offset_chan_file(1,tets)+2:offset_chan_file(filNum,tets)-offset_chan_file(1,tets)...
                        +naux2(filNum)+2)=...
                        aux_clu(:,offset_aux(filNum)+1+2:offset_aux(filNum)+naux2(filNum)+2);
                end
            else
                clu(:,:)=aux_clu(:,:);
            end
            
            %Limit to par.max_clus clusters
            for temper=1:szclu(1)
                too_large=find(clu(temper,3:end)>handles.par.max_clus-1);
                clu(temper,too_large+2)=-1;
            end
            
            %Define the classes
            class1=find(clu(temp,3:end)==0);
            class2=find(clu(temp,3:end)==1);
            class3=find(clu(temp,3:end)==2);
            class4=find(clu(temp,3:end)==3);
            class5=find(clu(temp,3:end)==4);
            try
                class0=setdiff(1:size(spikes,1), sort([class1; class2; class3; class4; class5]));
            catch
                class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
            end
            whos class*
            
            num_clusters = length(find([length(class1) length(class2) length(class3)...
                length(class4) length(class5) length(class0)] >= handles.par.min_clus));
            
            % If requested by user force membership
            %size(spikes,1)> handles.par.max_spk | ...
            if (handles.par.force_auto == 'y');
                classes = zeros(size(spikes,1),1);
                if length(class1)>=handles.par.min_clus; classes(class1) = 1; end
                if length(class2)>=handles.par.min_clus; classes(class2) = 2; end
                if length(class3)>=handles.par.min_clus; classes(class3) = 3; end
                if length(class4)>=handles.par.min_clus; classes(class4) = 4; end
                if length(class5)>=handles.par.min_clus; classes(class5) = 5; end
                f_in  = spikes(classes~=0,:);
                f_out = spikes(classes==0,:);
                class_in = classes(find(classes~=0),:);
                class_out = force_membership_wc(f_in, class_in, f_out, handles);
                classes(classes==0) = class_out;
                class0=find(classes==0);
                class1=find(classes==1);
                class2=find(classes==2);
                class3=find(classes==3);
                class4=find(classes==4);
                class5=find(classes==5);
            end
            
            %PLOTS
            figure(1)
            set(1,'Position',pos1);
            clf
            clus_pop = [];
            ylimit = [];
            subplot(2,5,6)
            temperature=handles.par.mintemp+temp*handles.par.tempstep;
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
            subplot(2,5,1)
            hold on
            drcluster=zeros(nspk,2);
            drcluster(:,2)= index';
            pooled_spikes=[spikes(:,:,1) spikes(:,:,2) spikes(:,:,3) spikes(:,:,4)];
            
            clus_pop = [clus_pop length(class0)];
            if length(class0) > 0;
                subplot(2,5,1);
                max_spikes=min(length(class0),handles.par.max_spikes);
                plot(pooled_spikes(class0(1:max_spikes),:)','k');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,5);
                hold on
                plot(pooled_spikes(class0(1:max_spikes),:)','k');
                plot(mean(pooled_spikes(class0,:),1),'c','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 0','Fontweight','bold')
                subplot(2,5,10)
                xa=diff(timestamp(class0));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                title([num2str(length(class0)) ' spikes']);
            end
            if length(class1) > handles.par.min_clus;
                clus_pop = [clus_pop length(class1)];
                subplot(2,5,1);
                max_spikes=min(length(class1),handles.par.max_spikes);
                plot(pooled_spikes(class1(1:max_spikes),:)','b');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,2);
                hold
                plot(pooled_spikes(class1(1:max_spikes),:)','b');
                plot(mean(pooled_spikes(class1,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 1','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,7)
                xa=diff(timestamp(class1));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','b','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                title([num2str(length(class1)) ' spikes']);
                drcluster(class1(:),1)=1;
            end
            
            if length(class2) > handles.par.min_clus;
                clus_pop = [clus_pop length(class2)];
                subplot(2,5,1);
                max_spikes=min(length(class2),handles.par.max_spikes);
                plot(pooled_spikes(class2(1:max_spikes),:)','r');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,3);
                hold
                plot(pooled_spikes(class2(1:max_spikes),:)','r');
                plot(mean(pooled_spikes(class2,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 2','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,8)
                xa=diff(timestamp(class2));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','r','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                drcluster(class2(:),1)=2;
                title([num2str(length(class2)) ' spikes']);
            end
            if length(class3) > handles.par.min_clus;
                clus_pop = [clus_pop length(class3)];
                subplot(2,5,1);
                max_spikes=min(length(class3),handles.par.max_spikes);
                plot(pooled_spikes(class3(1:max_spikes),:)','g');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,4);
                hold
                plot(pooled_spikes(class3(1:max_spikes),:)','g');
                plot(mean(pooled_spikes(class3,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 3','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,9)
                xa=diff(timestamp(class3));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','g','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                drcluster(class3(:),1)=3;
                title([num2str(length(class3)) ' spikes']);
            end
            if length(class4) > handles.par.min_clus;
                clus_pop = [clus_pop length(class4)];
                subplot(2,5,1);
                max_spikes=min(length(class4),handles.par.max_spikes);
                plot(pooled_spikes(class4(1:max_spikes),:)','c');
                xlim([1 size(pooled_spikes,2)]);
                drcluster(class4(:),1)=4;
            end
            if length(class5) > handles.par.min_clus;
                clus_pop = [clus_pop length(class5)];
                subplot(2,5,1);
                max_spikes=min(length(class5),handles.par.max_spikes);
                plot(pooled_spikes(class5(1:max_spikes),:)','m');
                xlim([1 size(pooled_spikes,2)]);
                drcluster(class5(:),1)=5;
            end
            
            % Rescale spike's axis
            if ~isempty(ylimit)
                ymin = min(ylimit(:,1));
                ymax = max(ylimit(:,2));
                
                if length(class1) > handles.par.min_clus; subplot(2,5,2); ylim([ymin ymax]); end
                if length(class2) > handles.par.min_clus; subplot(2,5,3); ylim([ymin ymax]); end
                if length(class3) > handles.par.min_clus; subplot(2,5,4); ylim([ymin ymax]); end
                if length(class0) > handles.par.min_clus; subplot(2,5,5); ylim([ymin ymax]); end
                
            end
            
            box off; hold on
            title([char(file_to_cluster) ' ch ' num2str(tets)],'Interpreter','none','Fontsize',14)
            features_name = handles.par.features;
            

            fig_name=[handles.directory '/' handles.drta_p.fullName(1:end-4) '/wave_fig2print_tetr' num2str(handles.drta_p.tets) handles.drta_p.fullName '.jpg'];
            saveas(figure(1),fig_name,'jpeg')
            
            cluster_class(timestamp_offset+1:timestamp_offset+szspikes(1),1:2) = drcluster(1:szspikes(1),1:2);
            szins=size(inspk);
            timestamp_offset=timestamp_offset+szspikes(1);
            outfile=['times_tetr' num2str(tets),'_' handles.drta_p.fullName(1:end-3) '.mat'];
            save(outfile, 'spikes','inspk', 'clu');
            
        end %if nspk>0
        %SAVE FILES
        cd(directory);
        
        %The joint_ file is used by wave_clus to process spikes using the output
        %produced here by cluster.exe
        outfile=['joint_' handles.drta_p.fullName(1:end-3) '.mat'];
        save(outfile, 'num_files', 'file_to_cluster','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
            'cluster_class','offset_for_chan', 'noSpikes');
        
        
        %Output data in jt_times files
        
        noSpikes=[];
        for filNum=1:handles.numfiles
            par = handles.par;
            
            
            handles.drta_p.draName=handles.files{filNum};
            handles.drta_p.fullName=handles.files{filNum};
            handles.par.filename=handles.files{filNum};
            
            
            
            load([handles.drta_p.draName(1:end-3),'.mat']);
            
            %If you want to pass drta_p variables upstream enter to drta_p
            drta_p.last_threshold=handles.drta_p.threshold;
            drta_p.f_strings=f_strings;
            
            save([handles.drta_p.draName(1:end-3),'.mat'],'data','params','drta_p');
            
            draq_p=params;
            draq_d=data;
            
            %The jt_times_ file is used by drs to analyze the data output by the
            %program
            outfile=['jt_times_' handles.drta_p.fullName(1:end-3) '.mat'];
            
            noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
            all_timestamp_per_file=[];
            cluster_class_per_file=[];
            for jj=1:4
                szat=size(all_timestamp_per_file);
                offset_for_chan(jj)=szat(2);
                all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
                cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
            end
            
            %If you want to pass drta_p variables upstream enter to drta_p
            handles.drta_p.last_threshold=handles.drta_p.threshold;
            drta_p.f_strings=f_strings;
            
            save(outfile, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d');
        end
    end %for tets=1:4
    
    
    
    %Output cluster results
    numclus=length(clus_pop)-1;
    
    
    outfileclus='cluster_results.txt';
    fout=fopen(outfileclus,'at+');
    fprintf(fout,'%s\t %s\t %g\t %d\t', file_to_cluster, features_name, temperature, numclus);
    for ii=1:numclus
        fprintf(fout,'%d\t',clus_pop(ii));
    end
    fprintf(fout,'%d\n',clus_pop(end));
    fclose(fout);
    
    
else
    h=msgbox('You must run drta before Do_wave_clusdrdg!!');
end

toc


% --- Executes on button press in load_run_rhd.
function load_run_rhd_Callback(hObject, eventdata, handles)
% hObject    handle to load_run_rhd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% %GET SPIKES prameters
tic

try
    close 1
catch
end

try
    close 2
catch
end
handles.file_name='';
handles.datatype ='CSC data'; %Needed to avoid output to figure
handles.flag=0;

print2file =1;                              %for saving printouts.
handles.par.prog = 'do_wave_clus';
handles.par=set_parameters_master(handles.par);
handles.par.dgorrhd=2;

directory=handles.directory;

bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');
pos1  = [bdwidth+(1/6)*scnsize(4),...
    (3/8)*scnsize(4),...
    scnsize(3)*(2/3) - 2*bdwidth,...
    scnsize(4)*(1/2) - (topbdwidth + bdwidth)];


set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8])

cd(directory);
num_files=handles.numfiles;


cd(directory);


handles.drta_p.draName=handles.files{1};
handles.drta_p.fullName=handles.files{1};
file_to_cluster = handles.files{1};
%     end
if ismac
    sub_directory=[directory '/' handles.drta_p.fullName(1:end-4)];
else
    sub_directory=[directory '\' handles.drta_p.fullName(1:end-4)];
end
load([handles.drta_p.draName(1:end-4),'.mat']);

an_error=0;
try
    handles.draq_p=params;
catch
    an_error=1;
end

if an_error==0
    handles.drta_p=drta_p;
    
    noSpikes=[];
    spikesExcluded=[];
    timestamp_offset=0;
    all_timestamp=[];   %Timestamps
    offset_for_chan=[]; %For each channel this has the offset within all_timestamp
    cluster_class=[];
    
    offset_chan_file=zeros(handles.numfiles,4);
    
    noSpikesChFl=zeros(handles.numfiles,4);
    noSpikesChFlExc=zeros(handles.numfiles,4);
    f_strings=cell(4,10);
    
    %Process data for each tetrode
    for tets=1:4
        tetrode_number=tets
        %For each tetrode find spikes, spike features and do processing
        
        spikes=[];
        index=[];
        timestamp=[];
        nexc=0;
        lbl='detect'
        
        
        
        for filNum=1:handles.numfiles
            file_number=filNum
            noSpikesChFlExc(filNum,tets)=0;
            handles.drta_p.draName=handles.files{filNum};
            cd(directory);
            load([handles.drta_p.draName(1:end-4),'.mat']);
            handles.draq_p=params;
            handles.par.sr=handles.draq_p.ActualRate;
            if handles.par.read_entire_file==1
                handles.draq_d=data;
                handles.draq_d.data=wc_read_dg(handles);
            end
            handles.drta_p=drta_p;
            handles.drta_p.draName=handles.files{filNum};
            handles.drta_p.fullName=handles.files{filNum};
            handles.par.filename=handles.files{filNum};
            handles.drta_p.tets=tets;
            sz_ts=size(timestamp);
            file_start_sz=sz_ts(2);
            
            for ii=1:handles.draq_d.noTrials
                %for ii=1:2
                trial_number=ii
                
                
                
                handles.drta_p.trialNo=ii;
                all_tet_ch_proc=1;
                for chno=4*(tets-1)+1:4*(tets-1)+4
                    
                    if (handles.drta_p.trial_ch_processed(chno,handles.drta_p.trialNo)==0)
                        all_tet_ch_proc=0;
                    end
                    
                end
                
                
                if (all_tet_ch_proc==1)
                    %Do not process if one of the channels for the tetrode was
                    %taken off for this trial
                    
                    [aux_spikes,thr,aux_index,aux_excluded,aux_nspk] = drdg_amp_detect_wc(handles);     % Detection with amp. thresh.
                    
                    if aux_nspk>0
                        shift=handles.draq_d.t_trial(ii)*handles.draq_p.ActualRate;
                        spikes=[spikes;aux_spikes];
                        index=[index aux_index];
                        timestamp=[timestamp (shift+aux_index)/handles.draq_p.ActualRate];
                        nexc=nexc+aux_excluded;
                        noSpikesChFlExc(filNum,tets)=noSpikesChFlExc(filNum,tets)+aux_excluded;
                    end
                end
                
                
            end
            sz_ts=size(timestamp);
            file_end_sz=sz_ts(2);
            offset_chan_file(filNum,tets)=timestamp_offset+file_start_sz;
            noSpikesChFl(filNum,tets)=file_end_sz-file_start_sz;
            %fclose(handles.drta_p.fullName);
        end
        
        
        
        nspk = size(spikes,1);
        %feature('memstats')
        all_timestamp(timestamp_offset+1:timestamp_offset+nspk)=timestamp(:);
        offset_for_chan(tets)=timestamp_offset;
        noSpikes(tets)=nspk;
        spikesExcluded(tets)=nexc;
        
        %Cluster starts here
        mkdir(sub_directory);
        cd(sub_directory);
        lbl='cluster'
        
        timestamp=timestamp*1000;
        handles.par.fname = ['data_' file_to_cluster];   %filename for interaction with SPC
        
        if nspk > handles.par.max_spk
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*handles.par.max_spk);
        else
            handles.par.min_clus = max(handles.par.min_clus_abs,handles.par.min_clus_rel*nspk);
        end
        
        % CALCULATES INPUTS TO THE CLUSTERING ALGORITHM.
        if nspk>0
            [inspk, f_str] = wave_features_drdg(spikes,handles);              %takes wavelet coefficients.
            for ii=1:10
                f_strings{tets,ii}=f_str{ii};
            end
            
            % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
            if size(spikes,1)> handles.par.max_spk
                naux = floor(handles.par.max_spk/handles.numfiles);
                inspk_aux=[];
                for filNum=1:handles.numfiles
                    if naux>noSpikesChFl(filNum,tets)
                        naux2(filNum)=noSpikesChFl(filNum,tets);
                    else
                        naux2(filNum)=naux;
                    end
                    szin=size(inspk_aux);
                    offset_aux(filNum)=szin(1);
                    inspk_aux = [inspk_aux; inspk(1+offset_chan_file(filNum,tets)-offset_chan_file(1,tets):...
                        offset_chan_file(filNum,tets)-offset_chan_file(1,tets)+naux2(filNum),:)];
                    
                end
            else
                inspk_aux = inspk;
            end
            
            
            %INTERACTION WITH SPC
            save(handles.par.fname_in,'inspk_aux','-ascii');
            [aux_clu, tree] = run_cluster(handles);
            [temp] = find_temp(tree,handles);
            
            
            %DEFINE CLUSTERS
            szclu=size(aux_clu);
            szspikes=size(spikes);
            clu=-1*ones(szclu(1),szspikes(1)+2);
            clu(:,1:2)=aux_clu(:,1:2);
            
            if size(spikes,1)> handles.par.max_spk
                for filNum=1:handles.numfiles
                    
                    clu(:,offset_chan_file(filNum,tets)+1-offset_chan_file(1,tets)+2:offset_chan_file(filNum,tets)-offset_chan_file(1,tets)...
                        +naux2(filNum)+2)=...
                        aux_clu(:,offset_aux(filNum)+1+2:offset_aux(filNum)+naux2(filNum)+2);
                end
            else
                clu(:,:)=aux_clu(:,:);
            end
            
            %Limit to par.max_clus clusters
            for temper=1:szclu(1)
                too_large=find(clu(temper,3:end)>handles.par.max_clus-1);
                clu(temper,too_large+2)=-1;
            end
            
            %Define the classes
            class1=find(clu(temp,3:end)==0);
            class2=find(clu(temp,3:end)==1);
            class3=find(clu(temp,3:end)==2);
            class4=find(clu(temp,3:end)==3);
            class5=find(clu(temp,3:end)==4);
            try
                class0=setdiff(1:size(spikes,1), sort([class1; class2; class3; class4; class5]));
            catch
                class0=setdiff(1:size(spikes,1), sort([class1 class2 class3 class4 class5]));
            end
            whos class*
            
            num_clusters = length(find([length(class1) length(class2) length(class3)...
                length(class4) length(class5) length(class0)] >= handles.par.min_clus));
            
            % If requested by user force membership
            %size(spikes,1)> handles.par.max_spk | ...
            if (handles.par.force_auto == 'y');
                classes = zeros(size(spikes,1),1);
                if length(class1)>=handles.par.min_clus; classes(class1) = 1; end
                if length(class2)>=handles.par.min_clus; classes(class2) = 2; end
                if length(class3)>=handles.par.min_clus; classes(class3) = 3; end
                if length(class4)>=handles.par.min_clus; classes(class4) = 4; end
                if length(class5)>=handles.par.min_clus; classes(class5) = 5; end
                f_in  = spikes(classes~=0,:);
                f_out = spikes(classes==0,:);
                class_in = classes(find(classes~=0),:);
                class_out = force_membership_wc(f_in, class_in, f_out, handles);
                classes(classes==0) = class_out;
                class0=find(classes==0);
                class1=find(classes==1);
                class2=find(classes==2);
                class3=find(classes==3);
                class4=find(classes==4);
                class5=find(classes==5);
            end
            
            %PLOTS
            figure(1)
            set(1,'Position',pos1);
            clf
            clus_pop = [];
            ylimit = [];
            subplot(2,5,6)
            temperature=handles.par.mintemp+temp*handles.par.tempstep;
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
            subplot(2,5,1)
            hold on
            drcluster=zeros(nspk,2);
            drcluster(:,2)= index';
            pooled_spikes=[spikes(:,:,1) spikes(:,:,2) spikes(:,:,3) spikes(:,:,4)];
            
            clus_pop = [clus_pop length(class0)];
            if length(class0) > 0;
                subplot(2,5,1);
                max_spikes=min(length(class0),handles.par.max_spikes);
                plot(pooled_spikes(class0(1:max_spikes),:)','k');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,5);
                hold on
                plot(pooled_spikes(class0(1:max_spikes),:)','k');
                plot(mean(pooled_spikes(class0,:),1),'c','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 0','Fontweight','bold')
                subplot(2,5,10)
                xa=diff(timestamp(class0));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                title([num2str(length(class0)) ' spikes']);
            end
            if length(class1) > handles.par.min_clus;
                clus_pop = [clus_pop length(class1)];
                subplot(2,5,1);
                max_spikes=min(length(class1),handles.par.max_spikes);
                plot(pooled_spikes(class1(1:max_spikes),:)','b');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,2);
                hold
                plot(pooled_spikes(class1(1:max_spikes),:)','b');
                plot(mean(pooled_spikes(class1,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 1','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,7)
                xa=diff(timestamp(class1));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','b','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                title([num2str(length(class1)) ' spikes']);
                drcluster(class1(:),1)=1;
            end
            
            if length(class2) > handles.par.min_clus;
                clus_pop = [clus_pop length(class2)];
                subplot(2,5,1);
                max_spikes=min(length(class2),handles.par.max_spikes);
                plot(pooled_spikes(class2(1:max_spikes),:)','r');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,3);
                hold
                plot(pooled_spikes(class2(1:max_spikes),:)','r');
                plot(mean(pooled_spikes(class2,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 2','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,8)
                xa=diff(timestamp(class2));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','r','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                drcluster(class2(:),1)=2;
                title([num2str(length(class2)) ' spikes']);
            end
            if length(class3) > handles.par.min_clus;
                clus_pop = [clus_pop length(class3)];
                subplot(2,5,1);
                max_spikes=min(length(class3),handles.par.max_spikes);
                plot(pooled_spikes(class3(1:max_spikes),:)','g');
                xlim([1 size(pooled_spikes,2)]);
                subplot(2,5,4);
                hold
                plot(pooled_spikes(class3(1:max_spikes),:)','g');
                plot(mean(pooled_spikes(class3,:),1),'k','linewidth',2)
                xlim([1 size(pooled_spikes,2)]);
                title('Cluster 3','Fontweight','bold')
                ylimit = [ylimit;ylim];
                subplot(2,5,9)
                xa=diff(timestamp(class3));
                [n,c]=hist(xa,0:1:100);
                bar(c(1:end-1),n(1:end-1))
                xlim([0 50])
                set(get(gca,'children'),'facecolor','g','linewidth',0.01)
                xlabel([num2str(sum(n(1:3))) ' in < 3ms'])
                drcluster(class3(:),1)=3;
                title([num2str(length(class3)) ' spikes']);
            end
            if length(class4) > handles.par.min_clus;
                clus_pop = [clus_pop length(class4)];
                subplot(2,5,1);
                max_spikes=min(length(class4),handles.par.max_spikes);
                plot(pooled_spikes(class4(1:max_spikes),:)','c');
                xlim([1 size(pooled_spikes,2)]);
                drcluster(class4(:),1)=4;
            end
            if length(class5) > handles.par.min_clus;
                clus_pop = [clus_pop length(class5)];
                subplot(2,5,1);
                max_spikes=min(length(class5),handles.par.max_spikes);
                plot(pooled_spikes(class5(1:max_spikes),:)','m');
                xlim([1 size(pooled_spikes,2)]);
                drcluster(class5(:),1)=5;
            end
            
            % Rescale spike's axis
            if ~isempty(ylimit)
                ymin = min(ylimit(:,1));
                ymax = max(ylimit(:,2));
                
                if length(class1) > handles.par.min_clus; subplot(2,5,2); ylim([ymin ymax]); end
                if length(class2) > handles.par.min_clus; subplot(2,5,3); ylim([ymin ymax]); end
                if length(class3) > handles.par.min_clus; subplot(2,5,4); ylim([ymin ymax]); end
                if length(class0) > handles.par.min_clus; subplot(2,5,5); ylim([ymin ymax]); end
                
            end
            
            box off; hold on
            title([char(file_to_cluster) ' ch ' num2str(tets)],'Interpreter','none','Fontsize',14)
            features_name = handles.par.features;
            
            fig_name=[handles.directory '/' handles.drta_p.fullName(1:end-4) '/wave_fig2print_tetr' num2str(handles.drta_p.tets) handles.drta_p.fullName '.jpg'];
            saveas(figure(1),fig_name,'jpeg')
            
            cluster_class(timestamp_offset+1:timestamp_offset+szspikes(1),1:2) = drcluster(1:szspikes(1),1:2);
            szins=size(inspk);
            timestamp_offset=timestamp_offset+szspikes(1);
            outfile=['times_tetr' num2str(tets),'_' handles.drta_p.fullName(1:end-4) '.mat'];
            save(outfile, 'spikes','inspk', 'clu');
            
        end %if nspk>0
        %SAVE FILES
        cd(directory);
        
        %The joint_ file is used by wave_clus to process spikes using the output
        %produced here by cluster.exe
        outfile=['joint_' handles.drta_p.fullName(1:end-4) '.mat'];
        save(outfile, 'num_files', 'file_to_cluster','directory','sub_directory','offset_chan_file', 'noSpikesChFl','spikesExcluded','all_timestamp',...
            'cluster_class','offset_for_chan', 'noSpikes');
        
        
        %Output data in jt_times files
        
        noSpikes=[];
        for filNum=1:handles.numfiles
            par = handles.par;
            
            
            handles.drta_p.draName=handles.files{filNum};
            handles.drta_p.fullName=handles.files{filNum};
            handles.par.filename=handles.files{filNum};
            
            
            
            load([handles.drta_p.draName(1:end-4),'.mat']);
            
            %If you want to pass drta_p variables upstream enter to drta_p
            drta_p.last_threshold=handles.drta_p.threshold;
            drta_p.f_strings=f_strings;
            
            save([handles.drta_p.draName(1:end-4),'mat'],'data','params','drta_p');
            
            draq_p=params;
            draq_d=data;
            
            %The jt_times_ file is used by drs to analyze the data output by the
            %program
            outfile=['jt_times_' handles.drta_p.fullName(1:end-4) '.mat'];
            
            noSpikes(1,1:4)=noSpikesChFl(filNum,1:4);
            all_timestamp_per_file=[];
            cluster_class_per_file=[];
            for jj=1:4
                szat=size(all_timestamp_per_file);
                offset_for_chan(jj)=szat(2);
                all_timestamp_per_file=[all_timestamp_per_file all_timestamp(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
                cluster_class_per_file=[cluster_class_per_file cluster_class(offset_chan_file(filNum,jj)+1:offset_chan_file(filNum,jj)+noSpikes(jj))];
            end
            
            %If you want to pass drta_p variables upstream enter to drta_p
            handles.drta_p.last_threshold=handles.drta_p.threshold;
            drta_p.f_strings=f_strings;
            
            save(outfile, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d');
        end
    end %for tets=1:4
    
    
    
    %Output cluster results
    numclus=length(clus_pop)-1;
    
    
    outfileclus='cluster_results.txt';
    fout=fopen(outfileclus,'at+');
    fprintf(fout,'%s\t %s\t %g\t %d\t', file_to_cluster, features_name, temperature, numclus);
    for ii=1:numclus
        fprintf(fout,'%d\t',clus_pop(ii));
    end
    fprintf(fout,'%d\n',clus_pop(end));
    fclose(fout);
    
    
else
    h=msgbox('You must run drta before Do_wave_clusdrdg!!');
end

toc



% --- Executes on button press in load_and_run.
function load_and_run_Callback(hObject, eventdata, handles)
% hObject    handle to load_and_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[files,path]=uigetfile({'*.dg';'*.rhd'}, 'Select dg or rhd files for wave_clus batch analysis', 'MultiSelect', 'on');
handles.directory=path(1:end-1);

if ischar(files)
    %This is one file
    handles.files{1}=files;
    handles.numfiles=1;
else
    handles.files=files;
    handles.numfiles=length(handles.files);
end

first_file=handles.files{1};
 
set(handles.whichFile,'String',first_file);

if strcmp(first_file(end-2:end),'rhd')
    load_run_rhd_Callback(hObject, eventdata, handles);
end

if strcmp(first_file(end-2:end),'.dg')
    load_run_dg_Callback(hObject, eventdata, handles);
end
