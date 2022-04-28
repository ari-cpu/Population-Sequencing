% input: _bcfilt.vcf files

%% !!!!!!!!!!!!! remember to change the name for saving !!!!!!!

clearvars

path ="/home/isabel/shared_Master/Matlab/"; %path where output summary is saved

filepath = "/home/isabel/MasterBackup/"; %path to input files

filesuffix="_bcfilt.vcf";

mode='ad'; % filter by AD thr or DP and variant frequency ?


% set thresholds
cutoff_dp=900;
cutoff_freq=0.02;
cutoff_ad=20; %!!!!!!!!!!!!!!!!!!!!!!!!!!anpassen je nach DP

sampleName=["LibSCBval89_bcfmp"] %input sample names
for i=1:length(sampleName)

    
%open bcfilt file and assigne cells 
fid=fopen(filepath+sampleName(i)+filesuffix);

imp=textscan(fid,'%*s %d . %s %s %*f . %*s %*s %s', 'CommentStyle', '#');

clear fid

 pos=imp(1);
 ref=imp(2);
 alt=imp(3);
 countinfo=imp(4);



% counts refers to the individual counts for every possible allele at a given
% position
countinfo=cellfun(@(a) split(a, ':'), countinfo, 'UniformOutput', false);
countinfo=countinfo{:};
counts=countinfo(:, 2);

counts=cellfun(@(a) split(a, ','), counts,'UniformOutput', false);



% format counts, calculate total count at each position (aka DP); seperate
% alt counts
counts=cellfun(@str2double, counts, 'UniformOutput', false);
sumcounts=cellfun(@sum, counts, 'UniformOutput', false);
altcounts=cellfun(@(a) a(2:end), counts, 'UniformOutput', false);
refcounts=cellfun(@(a) a(1), counts, 'UniformOutput', false);


%h_dp=histogram(cell2mat(sumcounts), 100);
%xlabel('Read Depth');
%legend("Read depth per base" + newline + "Mean: " + string(mean(cell2mat(sumcounts)))+...
%    newline + "Std: " + string(std(cell2mat(sumcounts))));
%% first filter:DP
mask_dp=cell2mat(sumcounts)>=cutoff_dp;


%% second filter AD
%find the maximum alternative AD at each position and filter 
[max, maxidx]=cellfun(@max, altcounts);
%ad_hist=histogram(max, 200);
%xlabel('AD of most abundant alternate');
%xlim([0, 50]);
%leg="Distribution of the AD of the most abundant alternate"; 
%mleg="Mean: " + string(mean(max));
%stdleg="std: " + string(std(max));
%legend(leg);




mask_ad=max>=cutoff_ad;

%calculate frequency for the maximum alternative allele at each position
freq=cellfun(@(a, b) a/b, num2cell(max), sumcounts, 'UniformOutput', false);

%h_freq=histogram(cell2mat(freq), 200);
%xlabel('variance frequency');
%xlim([0, 0.05]);
%legend({'Distribution of variance frequencies'});

mask_fr=cell2mat(freq)>=cutoff_freq;


%% apply filters depending on mode (ad or dp)

if mode=='dp'

mask_tot=mask_dp & mask_fr;

elseif mode=='ad';
    
    mask_tot=mask_ad;
end    

%apply the chosen filter to all other variables (reference allele count, alt allele count, position)
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
 
%fill Summary 
SNPSummary(i).FilterSummary.filtered=cell2struct(horzcat(num2cell([posfilt(:)]),refFilt, ...
    altfilt, freq(mask_tot), sumcounts(mask_tot), refcounts(mask_tot)),["pos", "ref", "alt", "varfreq","DP", "vfRef"], 2);




SNPSummary(i).Sample=sampleName(i);

clear mask_tot 

end

%save summmary
save(path+datestr(now, 'yyyymmdd')+"_LibSCBval89_BCFSummary.mat", 'SNPSummary');















