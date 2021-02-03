# drgWaveClus

drgWaveClus is a modification of the wave_clus pacakge of Quian-Quiroga to sort spikes using superparamagnetic clustering

Quiroga, R.Q., Nadasdy, Z., and Ben Shaul, Y. (2004). Unsupervised spike detection and sorting with wavelets and superparamagnetic clustering. Neural Comput 16, 1661-1687.
Quiroga, R.Q., and Panzeri, S. (2009). Extracting information from neuronal populations: information theory and decoding approaches. Nature Reviews Neuroscience 10, 173-185.


## Manuscripts

N/A

## Toolboxes

You should install the following MATLAB toolboxes:

Wavelet


## Spike sorting

Run Do_wave_clus with INTAN output files with metadata generated by running drta

Note: Neither the name of the directory nor the file name can have spaces because there will be an error when transferring the data to cluster.exe

Note: In MacOS Catalina I get the error: /bin/bash: ./cluster_maci.exe: Bad CPU type in executable

Then run wave_clus for each tetrode

## Output

The output is saved in a joint_name.mat file

offset_for_chan has the offset in all_timestamp for spikes for each tetrode

all_timestamp is a vector with all spike times in seconds

cluster_class is a two column vector. The first column has the ID for each spike. 0 were the rejected spikes, 1 to n are the sorted spikes. The second column has the index for the spike in the rhd file. Note that timestamp=index/acquisition rate.

 
