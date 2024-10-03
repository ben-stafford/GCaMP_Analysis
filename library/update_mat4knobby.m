%update_mat4knobby
%currently runs for all files in the current directory, specified by "dayx_xxx_xxx"
files=dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}))
currentFolder = pwd;

for i=1:length(directoryNames)
    str=sprintf('%s',directoryNames{i});
    cd(str)
    subfiles= dir; 
    fileNames= {subfiles.name};
    fileNames= fileNames(~ismember(fileNames,{'.','..'}));
    
    tf=contains('day',fileNames);

    mat_dir=dir('*day*.mat'); dot=find(mat_dir.name=='.'); fname=mat_dir.name(1:dot-1); %fname is the file name.
    load(fname);
    info.config.knobby.pos.x=[]; info.config.knobby.pos.y=[]; info.config.knobby.pos.z=[]; info.config.knobby.pos.a=[];
    save(mat_dir.name,'info')

    cd ..
end