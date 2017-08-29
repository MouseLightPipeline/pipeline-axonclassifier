sample = '2017-08-10'
h5folder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s/classifier_output',sample)
tiffolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',sample)
logfolder = '/groups/mousebrainmicro/mousebrainmicro/LOG/pipeline/'
axlogfiles = dir(fullfile(logfolder,['ax-2017-08*.txt']));
%%
clear args
tiffilelist = fullfile(tiffolder,'listtiffiles');
fid = fopen(tiffilelist);
inputlist = textscan(fid,'%s','Delimiter','\n');inputlist=inputlist{1};
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
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 10;
%%
% for every tif file, check if there is a h5 file and it is complete
numtif = length(inputlist);
myshfile = 'missing.sh'
fid = fopen(myshfile,'w');
for ii=26972:numtif
    inputname = inputlist{ii};
    [subpath,name,ext]=fileparts(inputname(length(tiffolder)+1:end));
    % check h5
    n1 = strrep(name,'ngc','prob');
    n2 = strrep(n1,'tif','h5');
    n3 = strrep(n1,'desc','txt');
    h5name = fullfile(h5folder,subpath,[n2,'.h5']);
    if exist(h5name,'file')
        % check completion
        spC = strsplit(subpath,'/');
        logfile = ['ax-',spC{2},'-',spC{3},'-',spC{4},n1(end-1:end),'.txt'];
        logpath = fullfile(logfolder,logfile);
        [q,w] = system(['grep ','"Completed Batch Processing"',' ',logpath]);
        if q
            q
%             % break
%         elseif 0
            %%
            % append to bsub
            randString = s( ceil(rand(1,sLength)*numRands) );
            jobname = sprintf('ilp_%05d-%s',ii,randString);
            
            infiles=inputname;
            % TODO: rename using prob
            subname = strsplit(name,{'-','.'});
            if ~exist(fullfile(h5folder,subpath), 'dir')
                mkdir(fullfile(h5folder,subpath));
                unix(sprintf('chmod g+rwx %s',fullfile(h5folder,subpath)));
            end
            outputformat=fullfile(h5folder,subpath,sprintf('%s-%s.%s.%s',subname{1},nametag,subname{3},outext));
            argsout = sprintf('''%s --headless  --cutout_subregion="[(None,None,None,0),(None,None,None,1)]" --logfile=%s --project=%s --output_format=%s --output_filename_format=%s %s''',...
                ilastikloc,logpath,ilpfile,outextformat,outputformat,infiles);
            % make sure name doesnot have any '.'
            name(name=='.')=[];
            %     mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d qsub -pe batch %d -l d_rt=1400 -N t-%d-%s -j y -o /dev/null -b y -cwd -V %s\n',numcores,memsize,numcores,ii,jobname,argsout);
            mysub = sprintf('LAZYFLOW_THREADS=%d LAZYFLOW_TOTAL_RAM_MB=%d bsub -n%d -We 25 -J t-%d-%s -o /dev/null %s\n',numcores,memsize,numcores,ii,jobname,argsout);
            fwrite(fid,mysub);
        end
    end
end
fclose(fid);

