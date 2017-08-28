% read target folder 
inputTifFolder = '/groups/mousebrainmicro/mousebrainmicro/data/2017-08-10/Tiling'

clear args
pathfile = fullfile(inputTifFolder,'listtiffiles');
args.ext = 'tif';
args.level = 3;
if exist(pathfile, 'file') == 2
    % load file directly
else
    args.fid = fopen(pathfile,'w');
    recdir(inputTifFolder,args)
    % make it write protected
    unix(sprintf('chmod -w %s',pathfile));
end
%%
filename = pathfile;
fid = fopen(filename);
targetlist = textscan(fid,'%s','Delimiter','\n');targetlist=targetlist{1};
fclose(fid);
%%
inputClassifierFolder = '/nrs/mouselight/cluster/classifierOutputs/2017-08-10/classifier_output'
clear args
pathfile = fullfile(inputClassifierFolder,'listh5files');
args.ext = 'h5';
args.level = 3;
if exist(pathfile, 'file') == 2
    % load file directly
else
    args.fid = fopen(pathfile,'w');
    recdir(inputClassifierFolder,args)
    % make it write protected
    unix(sprintf('chmod -w %s',pathfile));
end
%%

filename = pathfile;
fid = fopen(filename);
h5list = textscan(fid,'%s','Delimiter','\n');h5list=h5list{1};
fclose(fid);
%%

tiffolder = '/groups/mousebrainmicro/mousebrainmicro/data/2017-08-10/Tiling/';
h5folder = '/nrs/mouselight/cluster/classifierOutputs/2017-08-10/classifier_output/'
numTiles = length(targetlist);
keepit = zeros(length(h5list),1);
for ii=1:numTiles
    
    %%
    % check if target files are ready
    nametmp = targetlist{ii}(length(tiffolder):end);
    n1 = strrep(nametmp,'ngc','prob');
    n2 = strrep(n1,'tif','h5');
    n3 = strrep(n1,'desc','txt');
    
    inputname = fullfile(h5file,n2);
    idx = find(strcmp(h5list, inputname)); % single line engine
    keepit(idx) = 1;
end
%%
backupfolder = '/nrs/mouselight/cluster/classifierOutputs/2017-08-10/classifier_backup' 
mkdir(backupfolder)
for jj = find(~keepit)'
    foldname = fileparts(h5list{jj}(length(h5folder):end));
    mkdir(fullfile(backupfolder,foldname))
    unix(sprintf('mv %s %s',fullfile(h5folder,foldname),fullfile(backupfolder,foldname)))
end


















