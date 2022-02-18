clearvars

path = "/home/isabel/Documents/shared_Master/Matlab/";

filepath ="/home/isabel/Documents/shared_Master/Matlab/bcf/mock/";

filesuffix="_bcfilt.vcf";


sampleName="Bsub_MockReads_Mix_changeAll_bcfmp";

fid=fopen(filepath+sampleName+filesuffix);

imp=textscan(fid,'%*s %d . %s %s %*f . %*s %*s %s', 'CommentStyle', '#');


 pos=imp(1);
 ref=imp(2);
 alt=imp(3);
 countinfo=imp(4);




countinfo=cellfun(@(a) split(a, ':'), countinfo, 'UniformOutput', false);
countinfo=countinfo{:};
counts=countinfo(:, 2);




counts=cellfun(@(a) split(a, ','), counts,'UniformOutput', false);
cutoff_dp=1000;
cutoff_freq=0.01;


% format counts
counts=cellfun(@str2double, counts, 'UniformOutput', false);
sumcounts=cellfun(@sum, counts, 'UniformOutput', false);
altcounts=cellfun(@(a) a(2:end), counts, 'UniformOutput', false);
%refcounts=cellfun(@(a) a(1), counts, 'UniformOutput', false);


h_dp=histogram(cell2mat(sumcounts));
%% first filter:DP
mask_dp=cell2mat(sumcounts)>=cutoff_dp;

%% second filter AD
[max, maxidx]=cellfun(@max, altcounts);
freq=cellfun(@(a, b) a/b, num2cell(max), sumcounts, 'UniformOutput', false);
mask_ad=cell2mat(freq)>=cutoff_freq;


%% apply filters 

mask_tot=mask_dp & mask_ad;


%refFilt=refcounts(mask_tot);
ref=ref{:};
refFilt=ref(mask_tot);
refFilt=cellfun(@convertCharsToStrings, refFilt, 'UniformOutput', false);


alt=cellfun(@(a) split(a, ','), alt{:},'UniformOutput', false);
altmax=cellfun(@(a, b) a{b}, alt, num2cell(maxidx), 'UniformOutput', false);
altfilt=altmax(mask_tot);
altfilt=cellfun(@convertCharsToStrings, altfilt, 'UniformOutput', false);
pos=pos{:};
posfilt=pos(mask_tot);
refFilt=cellfun(@convertCharsToStrings, refFilt, 'UniformOutput', false);




%h2=histogram(varfreq);
%SNPSummary.FilterSummary.filtered=cell2struct(vertcat(num2cell([snp(:).pos]),num2cell([snp(:).ref]),...
 %   maxalts,num2cell(varfreq))', ["pos", "ref", "alt", "varfreq"], 2);
SNPSummary.FilterSummary.filtered=cell2struct(horzcat(num2cell([posfilt(:)]),refFilt, ...
    altfilt, freq(mask_tot)), ["pos", "ref", "alt", "varfreq"], 2);




SNPSummary.Sample=sampleName;

save(path+datestr(now, 'yyyymmdd')+"_BCFSummary.mat", 'SNPSummary');
















