function directory=dr_set_directory()

%sets the directory
%WARNING: do not enter spaces in the name or the pathway
%this (sometimes) leads to an error
%
%
% For example:
% directory='/Users/restrepd/Documents/Projects/Gould/T1040080114_wt_set 1'
% with a space between 'set' and '1' lead to the following error:
%
% lbl =
% 
% cluster
% 
% chmod: /Users/restrepd/Documents/Projects/Gould/T1040080114_wt_set: No such file or directory
% chmod: 1/T1040080114/cluster_maci.exe: No such file or directory
% /bin/bash: ./cluster_maci.exe: Permission denied
% Error using load
% Unable to read file 'times_tetr1_T1040080114.dg_01.lab'. No
% such file or directory.
%
% whereas 
% directory='/Users/restrepd/Documents/Projects/Gould/T1040080114_wt_set1'
% does work

%directory='/Users/restrepd/Documents/Projects/Gould/T1040080114_wt_set1';
%directory='/Users/restrepd/Documents/Projects/Daniel/matlabprocessdatabatchprogram/acetotheylben/nolaser/159863';
directory='/Users/restrepd/Documents/Projects/Justin/INTANtest';
%directory='/Users/restrepd/Documents/Projects/Ming/6203/051817';
%directory='/Users/restrepd/Documents/Projects/INTAN/462017'