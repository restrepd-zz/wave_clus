function Plot_all_features(handles)
% function Plot_all_features(handles)

try
    close 11
catch
end

try
    close 12
catch
end

start_class=1;
USER_DATA = get(handles.wave_clus_figure,'userdata');
par = USER_DATA{1};
inspk = USER_DATA{7};
classes = USER_DATA{6};

colors = ['k' 'b' 'r' 'g' 'c' 'm' 'y' 'b' 'r' 'g' 'c' 'm' 'y' 'b'];

figure(11)
annotation('textbox',[0.22, 0.97, 0, 0],'string','1:2')
annotation('textbox',[0.85, 0.97, 0, 0],'string','1:10')
annotation('textbox',[0.85, 0.2, 0, 0],'string','9:10')
nclasses = max(classes);
inputs = size(inspk,2);

for i=1:inputs
    for j=i+1:inputs
        subplot(inputs,inputs,(i-1)*inputs+j)
        hold on
        for k=start_class:nclasses
            %this_clu=find(clusters==clu)
            class_aux = find(classes==k);
            %max_spikes = min(par.max_spikes,length(class_aux));
            %plot(inspk(class_aux(1:max_spikes),i),inspk(class_aux(1:max_spikes),j),['.' colors(k)],'markersize',.5)
            plot(inspk(class_aux,i),...
                inspk(class_aux,j),['.' colors(k+1)],'MarkerSize',1);
            axis off
        end
        
%         figure(14)
%         hold on
%         for k=0:nclasses
%             %this_clu=find(clusters==clu)
%             class_aux = find(classes==k);
%             %max_spikes = min(par.max_spikes,length(class_aux));
%             %plot(inspk(class_aux(1:max_spikes),i),inspk(class_aux(1:max_spikes),j),['.' colors(k)],'markersize',.5)
%             plot(inspk(class_aux,i),...
%                 inspk(class_aux,j),['.' colors(k+1)],'MarkerSize',1);
%             axis off
%         end
%        
%         pffft=1;
%         try
%             close 14
%         catch
%         end
    end
end

[coefs,scores,variances,t2] = princomp(inspk);
figure(12)
start_class=1;
for k=start_class:nclasses
    class_aux = find(classes==k);
    max_spikes = min(par.max_spikes,length(class_aux));
    plot3(scores(class_aux(1:max_spikes),1),scores(class_aux(1:max_spikes),2),scores(class_aux(1:max_spikes),3),['.' colors(k+1)],'markersize',.5)
    hold on
    axis off
end

title('3D PCA for all features')


