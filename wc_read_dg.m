function data=wc_read_dg(handles)


%
% Reads the data from the recording computer
%
%

if handles.draq_p.dgordra==1
    %This is dra, this has not been updated, and will not work
    trialNo=handles.p.trialNo
    offset=handles.draq_p.no_spike_ch*2*(sum(handles.draq_d.samplesPerTrial(1:handles.p.trialNo))-handles.draq_d.samplesPerTrial(handles.p.trialNo));
    fseek(fid,offset,'bof');
    data_vec=fread(fid,handles.draq_p.no_spike_ch*handles.draq_d.samplesPerTrial(handles.p.trialNo),'uint16');
    szdv=size(data_vec);
    data=reshape(data_vec,szdv(1)/handles.draq_p.no_spike_ch,handles.draq_p.no_spike_ch);
else
    %This is dg
    fid=fopen(handles.drta_p.draName,'r');
    no_unit16_per_ch=floor(handles.draq_p.sec_per_trigger*handles.draq_p.ActualRate);
    data=fread(fid,no_unit16_per_ch*handles.draq_p.no_chans*handles.draq_d.noTrials,'uint16');
    fclose(fid)
end