# README #

axonal enhancement classifier. It aims to normalize axon appereance from various scopes and labeling strategies by utilizing shape and appereance features. 

### TO RUN ###

* run deployIlastikonbrain.m

### INPUTS ###
deployIlastikonbrain(brain,tag,logfolder)
brain: sample name
tag(OPTIONAL): experiment name, useful to run multiple experiments on the same brain, default = ''
logfolder(OPTIONAL): output log folder. Needed to search for succesfull completion, default : <output_tile_folder>/.log

### USAGE ### 
deployIlastikonbrain(brain,tag,logfolder) function will replicate the input structure of inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain); into 
outputtfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s%s/',brain,tag);

unless specified, log files are created in <outputfolder>/.log folder witn logtag = 'ax-%s-log.%s.txt'
