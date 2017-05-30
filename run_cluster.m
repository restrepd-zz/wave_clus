function [clu, tree] = run_cluster(handles)

% The while... try was added because cluster_maci.exe sometime has a problem
% (?Segmentation fault?). Because of this while... runs cluster again until
% it works. Thus for example:
%
% SPC v2.1 (14/2/2000) -- efficient KNN and mstree
% SPC v2.0.1 (1/12/2000) -- StopRunAtBreak option
% SPC v2.0 (1/11/1999) -- new param file + revised data structures
% -------------------------------------------------------------------
% ./cluster_maci.exe times_tetr2_T7110032814.run: Segmentation fault
%
% SPC v2.1 (14/2/2000) -- efficient KNN and mstree
% SPC v2.0.1 (1/12/2000) -- StopRunAtBreak option
% SPC v2.0 (1/11/1999) -- new param file + revised data structures
% -------------------------------------------------------------------
% AverageInteraction  0.091856
% ClustersReported  12
% CharDist 2.246740
%
% And so on?
%IMPORTANT: Do not name directories with spaces.
%e.g. do not use: 'Problem 10-20'
%in that case you will get the following error in an infinite loop:
% chmod: /Users/restrepd/Documents/Projects/Jorge/Problem: No such file or directory
% chmod: 10-20/amwt6_oct152015_2/cluster_maci.exe: No such file or directory
% /bin/bash: ./cluster_maci.exe: Permission denied
%
% in that example name the directory: 'Problem_10_20' (I am supersticious, I am using underline here)



% not_a_problem=0;
% no_problems=0;
% while not_a_problem==0
%
%     try
dim=handles.par.inputs;
if handles.par.dgorrhd==1
    fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:end-3)];
else
    fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:end-4)];
end
fname_in=handles.par.fname_in;
%set(handles.file_name,'string','Running SPC ...');

% DELETE PREVIOUS FILES
save([fname '.dg_01.lab'],'dim','-ASCII');         delete([fname '.dg_01.lab']);
save([fname '.dg_01'],'dim','-ASCII');              delete([fname '.dg_01']);

dat=load(fname_in);
n=length(dat);
fid=fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname_in);
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'Dimensions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(handles.par.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(handles.par.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(handles.par.tempstep));
fprintf(fid,'SWCycles: %s\n',num2str(handles.par.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(handles.par.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
if num2str(handles.par.randomseed) ~= 0
    fprintf(fid,'ForceRandomSeed: %s\n',num2str(handles.par.randomseed));
end
fclose(fid);

[str,maxsize,endian]=computer;
handles.par.system=str;

 
switch handles.par.system
    case {'PCWIN','PCWIN64'}
        if exist([pwd '\cluster.exe'])==0
            directory = which('cluster.exe');
            copyfile(directory,pwd);
        end
        dos(sprintf('cluster.exe %s.run',fname));
    case {'MAC'}
        if exist([pwd '/cluster_mac.exe'])==0
            directory = which('cluster_mac.exe');
            copyfile(directory,pwd);
        end
        run_mac = sprintf('./cluster_mac.exe %s.run',fname);
        unix(run_mac);
    case {'MACI','MACI64'}
        if exist([pwd '/cluster_maci.exe'])==0
            directory = which('cluster_maci.exe');
            copyfile(directory,pwd);
        end
        %run_maci = sprintf([pwd '/cluster_maci.exe %s.run'],fname);
        system(['chmod u+x ' pwd '/cluster_maci.exe']);
        run_maci = sprintf('./cluster_maci.exe %s.run',fname);
        unix(run_maci);
    otherwise  %(GLNX86, GLNXA64, GLNXI64 correspond to linux)
        if exist([pwd '/cluster_linux.exe'])==0
            directory = which('cluster_linux.exe');
            copyfile(directory,pwd);
        end
        run_linux = sprintf('./cluster_linux.exe %s.run',fname);
        unix(run_linux);
end


clu=load([fname '.dg_01.lab']);
tree=load([fname '.dg_01']);
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname_in);
%         not_a_problem=1;
%     catch
%         no_problems=no_problems+1;
%         problem='You have too few spikes'
%     end
% end
