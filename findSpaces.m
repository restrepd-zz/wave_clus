function there_is_a_space = findSpaces(handles)

there_is_a_space=0;

for filNum=1:handles.numfiles
    this_file_name=handles.files{filNum};
    for ii=1:length(this_file_name)
        if strcmp(this_file_name(ii),' ')
            there_is_a_space=1;
        end
    end
end

this_directory=handles.directory;
for ii=1:length(this_directory)
    if strcmp(this_directory(ii),' ')
        there_is_a_space=1;
    end
end