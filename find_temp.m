function [temp] = find_temp(tree,handles);
% Selects the temperature.

num_temp=handles.par.num_temp;
min_clus=handles.par.min_clus;

aux =diff(tree(:,5));   % Changes in the first cluster size
aux1=diff(tree(:,6));   % Changes in the second cluster size
aux2=diff(tree(:,7));   % Changes in the third cluster size
aux3=diff(tree(:,8));   % Changes in the third cluster size

temp = 1;         % Initial value

% for t=1:num_temp-1;
%     % Looks for changes in the cluster size of any cluster larger than min_clus.
%     if ( aux(t) > min_clus | aux1(t) > min_clus | aux2(t) > min_clus | aux3(t) >min_clus )    
%         temp=t+1;         
%     end
% end

aux1=tree(:,6);

temp2=find(aux1==max(aux1));
if length(temp2)>1
    temp3=temp2(1);
    temp2=[];
    temp2=temp3;
end

if (aux1(temp2)>min_clus)
   temp=temp2; 
end

%In case the second cluster is too small, then raise the temperature a little bit 
if (temp == 1 & tree(temp,6) < min_clus)
    temp = 2;
end    
   
