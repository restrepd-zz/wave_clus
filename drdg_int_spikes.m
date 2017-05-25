function [spikes1] = drdg_int_spikes(spikes,handles)
%Interpolates with cubic splines to improve alignment.

w_pre=handles.par.w_pre;
w_post=handles.par.w_post;
ls = w_pre + w_post;

int_factor = handles.par.int_factor;
nspk=size(spikes,1);

s=1:size(spikes,2);
ints=1/int_factor:1/int_factor:size(spikes,2);

intspikes=zeros(1,length(ints));
spikes1=zeros(nspk,ls,4);

for (ii=1:4)
    this_ch=ii+4*(handles.drta_p.tets-1);
    %thrmax(ii)=handles.drta_p.upper_limit(this_ch);
    if (handles.drta_p.threshold(this_ch)>0)
        detect='pos';
        %thr(ii)=handles.drta_p.threshold(this_ch);
    else
        detect='neg';
        %thr(ii)=-handles.drta_p.threshold(this_ch);
    end
end

switch detect
    case 'pos'
        for i=1:nspk
            for elec=1:4
                intspikes(:) = spline(s,spikes(i,:,elec),ints);
                [maxi iaux]=max((intspikes(w_pre*int_factor:w_pre*int_factor+8)));
                iaux = iaux + w_pre*int_factor -1;
                spikes1(i,w_pre:-1:1,elec) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
                spikes1(i,w_pre+1:ls,elec) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
            end
        end
    case 'neg'
        for i=1:nspk
            for elec=1:4
                intspikes(:) = spline(s,spikes(i,:,elec),ints);
                [maxi iaux]=min((intspikes(w_pre*int_factor:w_pre*int_factor+8)));
                iaux = iaux + w_pre*int_factor -1;
                spikes1(i,w_pre:-1:1,elec) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
                spikes1(i,w_pre+1:ls,elec) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
            end
        end
    case 'both'
        for i=1:nspk
            for elec=1:4
                intspikes(:) = spline(s,spikes(i,:,elec),ints);
                [maxi iaux]=max(abs(intspikes(w_pre*int_factor:w_pre*int_factor+8)));
                iaux = iaux + w_pre*int_factor -1;
                spikes1(i,w_pre:-1:1,elec) = intspikes(iaux:-int_factor:iaux-w_pre*int_factor+int_factor);
                spikes1(i,w_pre+1:ls,elec) = intspikes(iaux+int_factor:int_factor:iaux+w_post*int_factor);
            end
        end
end
