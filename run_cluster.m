function [clu, tree] = run_cluster(handles)

% Does superparamagnetic clustering using cluster_mac.exe et al
%
% 
%IMPORTANT: Do not name directories with spaces.
% e.g. do not use: 'Problem 10-20'
% in that case you will get the following error
% chmod: /Users/restrepd/Documents/Projects/Jorge/Problem: No such file or directory
% chmod: 10-20/amwt6_oct152015_2/cluster_maci.exe: No such file or directory
% /bin/bash: ./cluster_maci.exe: Permission denied
%
% in that example name the directory: 'Problem_10_20' 
% 
% If you get this error the name of the file is too long!
% error:
% at line 74 of 'param.c': too long
% Error using load
% Unable to read file
% 'times_tetr1_20180502_5118_spm_iso_mo_180502_100240.dg_01'. No such file
% or directory.
% 
% Error in run_cluster (line 102)
% tree=load([fname '.dg_01']);

dim=handles.par.inputs;
if handles.par.dgorrhd==1
    %dg
    if length(handles.par.filename)<=27
        fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:end-3)];
    else
        fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:28)];
    end
else
    %rhd
    if length(handles.par.filename)<=28
        fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:end-4)];
    else
        fname=['times_tetr' num2str(handles.drta_p.tets) '_' handles.par.filename(1:28)];
    end
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

% Strange, when clu=load([fname... is run there is sometimes the following
% error: ./cluster_maci.exe times_tetr1_170530_144920_6203.run: Segmentation fault
% If one runs it again it does not happen. For now I have a while loop

keep_clustering=1;
ii_clust=0;
while keep_clustering==1
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
    
    try
        clu=load([fname '.dg_01.lab']);
        keep_clustering=0;
    catch
        ii_clust=ii_clust+1;
    end
    if ii_clust>4
        keep_clustering=0;
    end
end
tree=load([fname '.dg_01']);
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname_in);

