function [inputtiflist,inputh5list] = listfolders(sample)
% match folders between tif and h5
h5folder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/classifier_output',sample);
tiffolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',sample);
% list tif files
clear args
pathfile = fullfile(tiffolder,'listtiffiles');
args.ext = 'tif';
args.level = 3;
if exist(pathfile, 'file') == 2
    % load file directly
else
    args.fid = fopen(pathfile,'w');
    recdir(tiffolder,args)
    % make it write protected
    unix(sprintf('chmod -w %s',pathfile));
end
clear args
tiffilelist = fullfile(tiffolder,'listtiffiles');
fid = fopen(tiffilelist);
inputtiflist = textscan(fid,'%s','Delimiter','\n');inputtiflist=inputtiflist{1};
fclose(fid);

% list h5 files
clear args
pathfile = fullfile(h5folder,'listh5files');
args.ext = 'h5';
args.level = 3;
if exist(pathfile, 'file') == 2
    % load file directly
else
    args.fid = fopen(pathfile,'w');
    recdir(h5folder,args)
    % make it write protected
    unix(sprintf('chmod -w %s',pathfile));
end
clear args
h5filelist = fullfile(h5folder,'listh5files');
fid = fopen(h5filelist);
inputh5list = textscan(fid,'%s','Delimiter','\n');inputh5list=inputh5list{1};
fclose(fid);

% make sure to not include any additional h5 files

end
