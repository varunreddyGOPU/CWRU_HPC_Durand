function getSpikesCount(data)
load(data);
spikesCountName = regexprep(data, '.mat', '_count.txt')
fid = fopen( spikesCountName, 'wt' );
fprintf(fid, '%d', length(spikes));
fclose(fid);