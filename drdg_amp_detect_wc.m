function [spikes,thr,index, excluded, nspk] = drdg_amp_detect_wc(handles);
% Detect spikes with amplitude thresholding. Uses median estimation.
% Detection is done with filters set by fmin_detect and fmax_detect. Spikes
% are stored for sorting using fmin_sort and fmax_sort. This trick can
% eliminate noise in the detection but keeps the spikes shapes for sorting.
index=[];
spikes=[];
excluded=0;
sr=handles.par.sr;
w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
% w_pre=floor(handles.draq_p.ActualRate*handles.drta_p.dt_pre_snip);
% w_post=floor(handles.draq_p.ActualRate*handles.drta_p.dt_post_snip);
%ref=handles.par.ref;
ref=w_post;

% stdmin = handles.par.stdmin;
% stdmax = handles.par.stdmax;
% fmin_detect = handles.par.detect_fmin;
% fmax_detect = handles.par.detect_fmax;
% fmin_sort = handles.par.sort_fmin;
% fmax_sort = handles.par.sort_fmax;

[data, data_sub, exc_ch] = getTraceDatadrdg(handles);

 
if (handles.drta_p.doSubtract==1)
    tetr=handles.drta_p.tets;
    for jj=1:4
        if handles.drta_p.subtractCh(4*(tetr-1)+jj)<=18
            if handles.drta_p.subtractCh(4*(tetr-1)+jj)<=16
                %Subtract one of the channels
                
                    data1(:,jj)=data(:,jj)-data_sub(:,jj);  
                
            else
                if handles.drta_p.subtractCh(4*(tetr-1)+jj)==17
                    %Subtract tetrode mean
                    data1(:,jj)=data(:,jj)-mean(data(:,1:4),2);
                else
                    %Subtract mean of all data
                    data1(:,jj)=data(:,jj)-mean(data,2);
                end
            end
        end
    end
end

% FILTER THE DATA between 500 and 5000 Hz
d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',500,'HalfPowerFrequency2',5000, ...
    'SampleRate',floor(handles.draq_p.ActualRate));

if (handles.drta_p.doSubtract==1)
    xf =filtfilt(d,data1);
else
    xf =filtfilt(d,data);
end

%Now setup thresholds
set(handles.file_name,'string','Detecting spikes ...');

for (ii=1:4)
    this_ch=ii+4*(handles.drta_p.tets-1);
    thrmax(ii)=handles.drta_p.upper_limit(this_ch);
    if (handles.drta_p.threshold(this_ch)>0)
        detect='pos';
        thr(ii)=handles.drta_p.threshold(this_ch);
    else
        detect='neg';
        thr(ii)=-handles.drta_p.threshold(this_ch);
    end
end


% LOCATE SPIKE TIMES
switch detect
    case 'pos'
        nspk = 0;
        for elect=1:4
            xaux = find(xf(w_pre+2:end-w_post-2,elect) > thr(elect)) +w_pre+1;
            xaux0 = 0;
            for ii=1:length(xaux)
                if xaux(ii) >= xaux0 + ref
                    [maxi iaux]=max((xf(xaux(ii):xaux(ii)+floor(ref/2)-1,elect)));    %introduces alignment
                    nspk = nspk + 1;
                    index(nspk) = iaux + xaux(ii) -1;
                    elect_no(nspk)=elect;
                    xaux0 = index(nspk);
                end
            end
        end
        
    case 'neg'
        nspk = 0;
        for elect=1:4
           
            xaux = find(xf(w_pre+2:end-w_post-2, elect) < -thr(elect)) +w_pre+1;
            
            xaux0 = 0;
            for ii=1:length(xaux)
                if xaux(ii) >= xaux0 + ref
                    [maxi iaux]=min((xf(xaux(ii):xaux(ii)+floor(ref/2)-1,elect)));    %introduces alignment
                    nspk = nspk + 1;
                    index(nspk) = iaux + xaux(ii) -1;
                    elect_no(nspk)=elect;
                    xaux0 = index(nspk);
                end
            end
        end
    case 'both'
        nspk = 0;
        for elect=1:4
            xaux = find(abs(xf(w_pre+2:end-w_post-2,elect)) > thr(elect)) +w_pre+1;
            xaux0 = 0;
            for ii=1:length(xaux)
                if xaux(ii) >= xaux0 + ref
                    [maxi iaux]=max(abs(xf(xaux(ii):xaux(ii)+floor(ref/2)-1,elect)));    %introduces alignment
                    nspk = nspk + 1;
                    index(nspk) = iaux + xaux(ii) -1;
                    elect_no(nspk)=elect;
                    xaux0 = index(nspk);
                end
            end
        end
end
if nspk>0
    % Sort the putative spikes
    
    sortArray(:,1)=index';
    sortArray(:,2)=elect_no';
    sortArray=sortrows(sortArray,1);
    index=sortArray(:,1)';
    elect_no=sortArray(:,2)';
    szxf=size(xf);
    
    %Get rid of spikes that user would like excluded (for things like licks)
    %excluding parameters (handles.drta_p.exc_sn_thr) are entered in drta
    
    %In Anan's channelrhodopsin program this was not entered
    if ~isfield(handles.drta_p,'exc_sn')
        handles.drta_p.exc_sn=0;
    end
    
    if handles.drta_p.exc_sn==1
        %xexc = filtfilt(b,a,exc_ch);
        xexc = filtfilt(d,exc_ch);
        aux=[];
        jj=1;
        for ii=1:length(index)
            if index(ii) <= (w_pre+2)
                %Exclude if the spike is within w_pre of the start
                aux(jj)=ii;
                jj=jj+1;
            else
                if index(ii)>=szxf(1)-(w_post+2)
                    %Exclude if the spike is within w_post of the end
                    aux(jj)=ii;
                    jj=jj+1;
                    
                else
                    if handles.drta_p.exc_sn_thr>0
                        %Now exclude if this is within the exclude spikes list
                        
                        if sum((xexc(index(ii)-w_pre-2:index(ii)+w_post+2))>handles.drta_p.exc_sn_thr)>0
                            aux(jj)=ii;
                            jj=jj+1;
                        end  
                       
                    else
                        if sum((xexc(index(ii)-w_pre-2:index(ii)+w_post+2))<handles.drta_p.exc_sn_thr)>0
                            aux(jj)=ii;
                            jj=jj+1;
                        end 
                    end
                end
            end
        end
        index(aux)=[];
        elect_no(aux)=[];
        excluded=length(aux);
        nspk=length(index);
    end
    
    
    %Work through all spikes that were chosen in two different channels
    %Keep the largest spike
    ls=w_pre+w_post;
    aux=[];
    jj=1;
    for ii=1:length(index)
        %Is the next spike within the extent of this one?
        if ii+1<length(index)
            if index(ii+1)-index(ii) <= ls
                if index(ii+1)+w_post-2<=size(xf,1)
                    %Choose the larger spike
                    deltaii=max(xf(index(ii)-w_pre+2:index(ii)+w_post-2,elect_no(ii)))-min(xf(index(ii)-w_pre+2:index(ii)+w_post-2,elect_no(ii)));
                    deltaii1=max(xf(index(ii+1)-w_pre+2:index(ii+1)+w_post-2,elect_no(ii+1)))-min(xf(index(ii+1)-w_pre+2:index(ii+1)+w_post-2,elect_no(ii+1)));
                    if deltaii>deltaii1
                        aux(jj)=ii+1;
                        lastii=ii;
                        jj=jj+1;
                    else
                        aux(jj)=ii;
                        lastii=ii+1;
                        jj=jj+1;
                    end
                    no_steps=1;
                    
                    %See if the next one is the same event
                    if ii+2<length(index)
                        if index(ii+2)-index(ii) <= ls
                            if index(ii+2)+w_post-2<=size(xf,1)
                                deltaii=max(xf(index(lastii)-w_pre+2:index(lastii)+w_post-2,elect_no(lastii)))-min(xf(index(lastii)-w_pre+2:index(lastii)+w_post-2,elect_no(lastii)));
                                deltaii1=max(xf(index(ii+2)-w_pre+2:index(ii+2)+w_post-2,elect_no(ii+2)))-min(xf(index(ii+2)-w_pre+2:index(ii+2)+w_post-2,elect_no(ii+2)));
                                if deltaii>deltaii1
                                    aux(jj)=ii+2;
                                    jj=jj+1;
                                else
                                    aux(jj)=lastii;
                                    lastii=ii+2;
                                    jj=jj+1;
                                end
                                no_steps=2;
                                
                                %One more time see if the next one is the same event
                                if ii+3<length(index)
                                    if index(ii+3)-index(ii) <= ls
                                        if index(ii+3)+w_post-2<=size(xf,1)
                                            deltaii=max(xf(index(lastii)-w_pre+2:index(lastii)+w_post-2,elect_no(lastii)))-min(xf(index(lastii)-w_pre+2:index(lastii)+w_post-2,elect_no(lastii)));
                                            deltaii1=max(xf(index(ii+3)-w_pre+2:index(ii+3)+w_post-2,elect_no(ii+3)))-min(xf(index(ii+3)-w_pre+2:index(ii+3)+w_post-2,elect_no(ii+3)));
                                            if deltaii>deltaii1
                                                aux(jj)=ii+3;
                                                jj=jj+1;
                                            else
                                                aux(jj)=lastii;
                                                lastii=ii+3;
                                                jj=jj+1;
                                            end
                                            no_steps=3;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    
                    ii=ii+no_steps;
                end
            end
        end
    end
    index(aux)=[];
    elect_no(aux)=[];
    nspk=length(index);
    
    
    %Get rid of artifacts >=thrmax
    spikes=zeros(nspk,ls+4,4);
    xf=[xf; zeros(w_post,4)];
    
    jj=1;
    for i=1:nspk                          %Eliminates artifacts
        save_spike=1;
        for elec=1:4
            if max(abs( xf(index(i)-w_pre:index(i)+w_post,elec) )) >= thrmax(elec)
                save_spike=0;
            end
        end
        if save_spike==1
            spikes(jj,:,:)=xf(index(i)-w_pre-1:index(i)+w_post+2,:);
            jj=jj+1;
        end
    end
    if jj-1<nspk
        spikes(jj:nspk,:,:)=[];
        index(jj:nspk)=[];
    end
    
    % SPIKE STORING (with or without interpolation)
    switch handles.par.interpolation
        case 'n'
            spikes(:,end-1:end)=[];       %eliminates borders that were introduced for interpolation
            spikes(:,1:2)=[];
        case 'y'
            %Does interpolation
            spikes = drdg_int_spikes(spikes,handles);
    end
    
    if ~(strcmp(handles.datatype,'CSC data') & strcmp(handles.par.tmax,'all'))
        USER_DATA = get(handles.wave_clus_figure,'userdata');
        USER_DATA{2}=spikes;
        USER_DATA{3}=index*1000/sr;
        set(handles.wave_clus_figure,'userdata',USER_DATA);
        Plot_continuous_data(xf,handles,thr,thrmax)
    elseif handles.flag == 1
        Plot_continuous_data(xf(1:floor(60*sr)),handles,thr,thrmax)
    end
end
