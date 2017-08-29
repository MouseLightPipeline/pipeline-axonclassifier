%% check missing files after classifier
logfolder = '/groups/mousebrainmicro/mousebrainmicro/LOG/pipeline/'
inputTifFolder = '/groups/mousebrainmicro/mousebrainmicro/data/2017-08-10/Tiling'
axlogfiles = dir(fullfile(logfolder,['ax-2017-08*.txt']));

%%
nametag = 'prob'
numcores = 4;
memsize = numcores*7.5*1000;
ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.9-Linux/run_ilastik.sh'
ilpfile = '/groups/mousebrainmicro/mousebrainmicro/erhan_dm11/AxonClassifier/axon_uint16.ilp';
if 0
    outext = 'tif'
    outextformat = '"multipage tiff"'
else
    outext = 'h5'
    outextformat = '"hdf5"'
end
%%
logfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/LOG/%s%s/',brain,'offline');
mkdir(logfolder)

s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 10;

numTiles = length(targetlist);
% chech number of h5 files
myshfile = 'missing.sh'
fid = fopen(myshfile,'w');


parfor ii=1:length(axlogfiles)
    myfile = fullfile(logfolder,axlogfiles(ii).name);
    [~,mf,~] = fileparts(myfile);
    [q,w] = system(['grep ','"Completed Batch Processing"',' ',myfile]);
    finishedaxfiles(ii)=~q;
    if q
        % check if tif folder exists
        spT = strsplit(mf,'-');
        targetfold = fullfile(inputTifFolder,[spT{2},'-',spT{3},'-',spT{4}],spT{5},spT{6}(1:end-2));
        if exist(fullfile(targetfold,[spT{6},'.tif']),'file')
            % missing
            
            randString = s( ceil(rand(1,sLength)*numRands) );
            jobname = sprintf('ilp_%05d-%s',ii,randString);
            logfile=fullfile(logfolder,sprintf('ilp_%05d-%s.txt',ii,randString));
            infiles=targetlist{ii};
            
            [subpath,name,ext]=fileparts(targetlist{ii}(length(inputfolder)+1:end));
            % TODO: rename using prob
            subname = strsplit(name,{'-','.'});
            if ~exist(fullfile(out,subpath), 'dir')
                mkdir(fullfile(out,subpath));
                unix(sprintf('chmod g+rwx %s',fullfile(out,subpath)));
            end
            outputformat=fullfile(out,subpath,sprintf('%s-%s.%s.%s',subname{1},nametag,subname{3},outext));
            argsout = sprintf('''%s --headless  --cutout_subregion="[(None,None,None,0),(None,None,None,1)]" --logfile=%s --project=%s --output_format=%s --output_filename_format=%s %s''',...
                ilastikloc,logfile,ilpfile,outextformat,outputformat,infiles);
            % make sure name doesnot have any '.'
            name(name=='.')=[];
            %     mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l d_rt=1400 -N t-%d-%s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,ii,jobname,argsout);
            mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d bsub -n%d -We 25 -J t-%d-%s -o /dev/null %s\n',numcores,memsize,numcores,ii,jobname,argsout);
            fwrite(fid,mysub);
        end
    end
    
end

%%
% read target folder 

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
%%
filename = '/nrs/mouselight/cluster/classifierOutputs/2017-08-10/classifier_output/listh5files';
fid = fopen(filename);
h5list = textscan(fid,'%s','Delimiter','\n');h5list=h5list{1};
fclose(fid);
%%
nametag = 'prob'
numcores = 4;
memsize = numcores*7.5*1000;
ilastikloc = '/groups/mousebrainmicro/mousebrainmicro/cluster/software/ilastik-1.1.9-Linux/run_ilastik.sh'
ilpfile = '/groups/mousebrainmicro/mousebrainmicro/erhan_dm11/AxonClassifier/axon_uint16.ilp';
if 0
    outext = 'tif'
    outextformat = '"multipage tiff"'
else
    outext = 'h5'
    outextformat = '"hdf5"'
end
%%
logfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/LOG/%s%s/',brain,'offline');
mkdir(logfolder)

s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 10;

numTiles = length(targetlist);
% chech number of h5 files
myshfile = 'missing.sh'
fid = fopen(myshfile,'w');
for ii=9094:numTiles
    nametmp = targetlist{ii}(length(inputTifFolder)+1:end);
    n1 = strrep(nametmp,'ngc','prob');
    n2 = strrep(n1,'tif','h5');
    
    % check if classifier exists
    if exist(fullfile(h5folder,n2),'file')
        continue
    end
    
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    jobname = sprintf('ilp_%05d-%s',ii,randString);
    logfile=fullfile(logfolder,sprintf('ilp_%05d-%s.txt',ii,randString));
    infiles=targetlist{ii};
    
    [subpath,name,ext]=fileparts(targetlist{ii}(length(inputfolder)+1:end));
    % TODO: rename using prob 
    subname = strsplit(name,{'-','.'});
    if ~exist(fullfile(out,subpath), 'dir')
        mkdir(fullfile(out,subpath));
        unix(sprintf('chmod g+rwx %s',fullfile(out,subpath)));
    end
    outputformat=fullfile(out,subpath,sprintf('%s-%s.%s.%s',subname{1},nametag,subname{3},outext));
    argsout = sprintf('''%s --headless  --cutout_subregion="[(None,None,None,0),(None,None,None,1)]" --logfile=%s --project=%s --output_format=%s --output_filename_format=%s %s''',...
        ilastikloc,logfile,ilpfile,outextformat,outputformat,infiles);
    % make sure name doesnot have any '.'
    name(name=='.')=[];
%     mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l d_rt=1400 -N t-%d-%s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,ii,jobname,argsout);
    mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d bsub -n%d -We 25 -J t-%d-%s -o /dev/null %s\n',numcores,memsize,numcores,ii,jobname,argsout);
    fwrite(fid,mysub);    
    
%     n3 = strrep(n1,'desc','txt');

end






























