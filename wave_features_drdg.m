function [inspk,f_strings] = wave_features_drdg(aux_spikes,handles)
%Calculates the spike features

feature_points=handles.par.feature_points;
spikes=aux_spikes(:,1:feature_points,:);
scales = handles.par.scales;
feature = handles.par.features;
inputs = handles.par.inputs;
nspk=size(spikes,1);
ls = size(spikes,2);
%set(handles.file_name,'string','Calculating spike features ...');


% CALCULATES FEATURES
switch feature
    case {'wav','wavfpv'}
        cc=zeros(nspk,ls*4);
        for i=1:nspk
            for tetr=1:4
                % Wavelet decomposition before adding the if exist
                %             [c,l]=wavedec(spikes(i,:,tetr),scales,'haar');
                %             cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls);
                %After adding the if exist
                if exist('wavedec')                             % Looks for Wavelets Toolbox
                    % Wavelet decomposition
                    [c,l]=wavedec(spikes(i,:,tetr),scales,'haar');
                    cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls);
                    
                else
                    % Replaces Wavelets Toolbox, if not available
                    [c,l]=fix_wavedec(spikes(i,:,tetr),scales);
                    cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls);
                    
                end
            end
        end
%         for i=1:ls*4                                  % KS test for coefficient selection   
%             thr_dist = std(cc(:,i)) * 3;
%             thr_dist_min = mean(cc(:,i)) - thr_dist;
%             thr_dist_max = mean(cc(:,i)) + thr_dist;
%             aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i); %discard outliers
%             if length(aux) > 10;
%                 [ksstat]=test_ks(aux);
%                 sd(i)=ksstat;
%             else
%                 sd(i)=0;
%             end
%         end
%         [maxim ind]=sort(sd);
%         coeff(1:inputs)=ind(ls:-1:ls-inputs+1);
        
    case 'pca'
        cc=[];
        for tetr=1:4
            [C,S,L] = princomp(spikes(:,:,tetr));
            cc = [cc S];
            %         inputs = 3;
            %         coeff(1:3)=[1 2 3];
        end
        

        
    case 'all'
        cc=zeros(nspk,ls*4);
        for i=1:nspk                                % Wavelet decomposition
            for tetr=1:4
                %Before adding the if exist
               %[c,l]=wavedec(spikes(i,:,tetr),scales,'haar');
               %cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls);
                               
               %After adding the if exist
                if exist('wavedec')                             % Looks for Wavelets Toolbox
                    % Wavelet decomposition
                    [c,l]=wavedec(spikes(i,:,tetr),scales,'haar');
                    cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls);
                    
                else
                    % Replaces Wavelets Toolbox, if not available
                    [c,l]=fix_wavedec(spikes(i,:,tetr),scales);
                    cc(i,1+ls*(tetr-1):ls*tetr)=c(1:ls); 
                end
               
           end
       end

        %Now add PCA
        for tetr=1:4
            [C,S,L] = princomp(spikes(:,:,tetr));
            cc(:,ls*4+1+3*(tetr-1):ls*4+3+3*(tetr-1)) = S(:,1:3);
        end
        
        %Then add peak to valley
        for tetr=1:4
            cc(:,4*(ls+3)+tetr)=max(spikes(:,1:feature_points,tetr)')-min(spikes(:,1:feature_points,tetr)');
        end


end   

        % KS test for coefficient selection
        sizecc=size(cc);
        for i=1:sizecc(2)                                     
            thr_dist = std(cc(:,i)) * 3;
            thr_dist_min = mean(cc(:,i)) - thr_dist;
            thr_dist_max = mean(cc(:,i)) + thr_dist;
            aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i); %discard outliers
            if length(aux) > 10;
                [ksstat]=test_ks(aux);
                sd(i)=ksstat;
            else
                sd(i)=0;
            end
        end
        [maxim ind]=sort(sd);
        coeff(1:inputs)=ind(sizecc(2):-1:sizecc(2)-inputs+1);
        
        
        
        %plot histograms of the 10 coefficients
        figure(2)
        for jj=1:inputs
            subplot(ceil(sqrt(inputs)),ceil(sqrt(inputs)),jj);
            hist(cc(:,coeff(jj)),20);
          
            if coeff(jj)<=ls*4
                %This is a wavelength
                wavno=coeff(jj)-ls*floor(coeff(jj)/ls);
                tetno=ceil(coeff(jj)/ls);
                title(['w' num2str(wavno) ' e' num2str(tetno)],'Color','r');
                f_strings{jj}=['w' num2str(wavno) ' e' num2str(tetno)];
            end
               
            if (coeff(jj)>ls*4)&(coeff(jj)<(ls*4+3*4))
                %This is a PCA
                tetno=floor((coeff(jj)-4*ls)/3)+1;
                pcano=coeff(jj)-4*ls-(tetno-1)*3;
                title(['p' num2str(pcano) ' e',num2str(tetno)],'Color','r');
                f_strings{jj}=['p' num2str(pcano) ' e',num2str(tetno)];
            end
            
            if (coeff(jj)>(ls*4+3*4))
               %This is the PV
               tetno=coeff(jj)-4*ls-3*4;
               title(['p-v e',num2str(tetno) ],'Color','r');
               f_strings{jj}=['p-v e',num2str(tetno) ];
            end


        end
        set(gcf,'papertype','usletter','paperorientation','portrait','paperunits','inches')
        set(gcf,'paperposition',[.25 .25 10.5 7.8])
%         directory=dr_set_directory();
%         cd(directory);
       
        eval(['print -djpeg40 feat_fig2print_tetr' num2str(handles.drta_p.tets) handles.drta_p.fullName '.jpg']);


%CREATES INPUT MATRIX FOR SPC
inspk=zeros(nspk,inputs);
for ii=1:nspk
    for jj=1:inputs
        inspk(ii,jj)=cc(ii,coeff(jj));
    end
end
% 



switch feature
    case 'wavfpv'
       inspk(1:nspk,inputs)=max(spikes(:,1:feature_points)')-min(spikes(:,1:feature_points)'); 
end

