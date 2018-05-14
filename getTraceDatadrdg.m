function [data, data_sub, exc_ch]=getTraceDatadrdg(handles)
% Gets data for one trial

if handles.par.dgorrhd==1
    %This is a dg file
    trialNo=handles.drta_p.trialNo;
    if handles.par.read_entire_file==1
        
        
        jj=1;
        
        data_this_trial=handles.draq_d.data(floor(handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger*handles.draq_p.no_chans*(trialNo-1)+1):...
            floor(handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger*handles.draq_p.no_chans*trialNo)-2000);
        
        data=[];
        for ii=1+4*(handles.drta_p.tets-1):4+4*(handles.drta_p.tets-1)
            data(:,jj)=data_this_trial(floor((ii-1)*handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)+1:...
                floor((ii-1)*handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)...
                +floor(handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)-2000);
            
            jj=jj+1;
        end
        
    else
        fid=fopen(handles.drta_p.draName,'rb');
        %Note: two bytes per sample (uint16)
        
        %This is dg
        
        bytes_per_native=2;     %Note: Native is unit16
        size_per_ch_bytes=handles.draq_p.sec_per_trigger*handles.draq_p.ActualRate*bytes_per_native;
        no_unit16_per_ch=size_per_ch_bytes/bytes_per_native;
        trial_offset=handles.draq_p.no_chans*size_per_ch_bytes*(trialNo-1);
        data_this_trial=[];
        for ii=1:handles.draq_p.no_chans
            fseek(fid, (ii-1)*size_per_ch_bytes+trial_offset, 'bof');
            data_this_trial(:,ii)=fread(fid,no_unit16_per_ch,'uint16');
        end
        
        fclose(fid);
        
        jj=1;
        data=[];
        for ii=1+4*(handles.drta_p.tets-1):4+4*(handles.drta_p.tets-1)
            data(:,jj)=data_this_trial(:,ii);
            data_sub(:,jj)=data_this_trial(:,handles.drta_p.subtractCh(ii));
            jj=jj+1;
        end
    end
    
    
    
    exc_ch=[];
    
    %Now if there is an exclusion channel get the data from that channel
    
    %In Anan's channelrhodopsin program this was not entered
    if ~isfield(handles.drta_p,'exc_sn')
        handles.drta_p.exc_sn=0;
    end
    
    if handles.drta_p.exc_sn==1
        
        
        %This is now using the lick data to exclude lick artifacts
        if handles.par.read_entire_file==1
            ii=19; %lick channel
            lick_ch=data_this_trial(floor((ii-1)*handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)+1:...
                floor((ii-1)*handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)...
                +floor(handles.draq_p.ActualRate*handles.draq_p.sec_per_trigger)-2000);
        else
            lick_ch=data(:,19);
        end
        exc_ch=abs(lick_ch(1:end-1)-lick_ch(2:end));
        
        
        exc_ch=[exc_ch;0];
    end
    
else
    %rhd file
    data_this_trial=wcGetTraceDataRHD(handles);
    jj=1;
    data=[];
    for ii=1+4*(handles.drta_p.tets-1):4+4*(handles.drta_p.tets-1)
        data(:,jj)=data_this_trial(:,ii);
        data_sub(:,jj)=data_this_trial(:,handles.drta_p.subtractCh(ii));
        jj=jj+1;
    end
    
    exc_ch=[];
    
    %Now if there is an exclusion channel get the data from that channel
    
    %In Anan's channelrhodopsin program this was not entered
    if ~isfield(handles.drta_p,'exc_sn')
        handles.drta_p.exc_sn=0;
    end
    
    if handles.drta_p.exc_sn==1
        lick_ch=data_this_trial(:,19);
        exc_ch=abs(lick_ch(1:end-1)-lick_ch(2:end));
        exc_ch=[exc_ch;0];
    end
    
end



