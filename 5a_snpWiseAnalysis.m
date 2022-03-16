%input: SNPSummary or BCFSummary, SNPSummary of monoclonal data 
%generates a CallSummary containing found SNPs compared to expected SNPs,
%false positives etc
%basis is the comparison with expected SNPs from monoclonal sequencing data
%of all variants in the mix


clearvars

                     %                      %
                    % %    Define input    % %
                     %                      %

donor = "Bval" 
%sample = "Bsub_MockReads_Mix_changeAll_bcfmp" 
saveSum=false;


SNPPath = "/home/isabel/shared_Master/Matlab/"; %path where SNP/BCFSummary is found
artefacts=SNPPath+"Artefacts_Bs1662BsNCe.vcf";

SNPName{1} = "20220314_Mix10_BCFSummary.mat"; %name of the Summary you want to analyze
%SNPName{1}="20220310_Mix1_BCFSummary.mat";
%SNPName{2} = "20211124_comp_SNPSummary.mat";
SNPName{2} = "20220314_newMixes_comp_SNPSummary.mat"; %name of the SNPSummary containing the monoclonal data for comparison
%SNPName{2} = "20220217_comp_SNPSummary.mat"; %with mutations


PopData=load(SNPPath+SNPName{1});
MonoData=load(SNPPath+SNPName{2});

PopData=PopData.SNPSummary;
MonoData=MonoData.SNPSummary;

%name_idx=[PopData.Sample]==sample;
i=size(PopData);

for name_idx=1:i(2)% loop over all SNP/BCFSummaries

SampleData=PopData(name_idx);

snpsPop=[SampleData.FilterSummary.filtered];


% Read in artefacts.vcf (Ancestor reads on dictionary)
if ~isempty(artefacts)
fid = fopen(artefacts);
imp = textscan(fid, '%s %f %s %s %s %f %s %s %s %s', 'HeaderLines', 33);
artefact.pos = imp{2};
artefact.alt = imp{5};
fclose(fid);
clear imp;
end 
%remove artefacts
for a = 1 : length(artefact.pos)
    idx = find([snpsPop(:).pos] == artefact.pos(a));
    if ~isempty(idx) %&& strcmp(artefact.alt(a), SampleData.FilterSummary.filtered(idx).alt)
        snpsPop(idx)=[];
       %snpsPopPos(idx) = [];
       %snpsPopFreq(idx) = [];
    end   
end 
clear idx a

%initialize summary 
Summary(length(MonoData))=struct();

%% comparison with the Monoclonal Data, loop over all samples in the monoclonal data 
for i=1:length(MonoData)
Summary(i).MonoclonalSample=MonoData(i).Sample;
snpsMono=MonoData(i).FilterSummary.filtered;

% remove artefacts from monoclonal Data 
for a = 1 : length(artefact.pos)
    idx = find([snpsMono(:).pos] == artefact.pos(a));
    if ~isempty(idx) %&& strcmp(artefact.alt(a), MonoData(i).FilterSummary.filtered(idx).alt)
       snpsMono(idx) = [];
    end   
end 
clear idx a

Summary(i).numberOfSNPs=length([snpsMono(:).pos]);

%check if positions in monoclonal data are found in population data
[where, there]=ismember([snpsMono(:).pos], [snpsPop(:).pos]);
there_idx=there(where);



%check for faulty snps (right position bur wrong allele)
checkSNPs=[snpsMono(where).alt]==[snpsPop(there_idx).alt];

if  isempty(nonzeros(checkSNPs))
    Summary(i).NoFalseSNPs=0;
    Summary(i).FalseSNPS=[];
else
    Summary(i).NoFalseSNPs=length(nonzeros(checkSNPs(~(nonzeros(checkSNPs)))));
    Summary(i).FalseSNPs=snpsPop(~checkSNPs);
end

%compute ratio of snps that is found in population data
Summary(i).RatioFound=length(nonzeros(where))/length([snpsMono(:).pos]);



  

%list all found SNPs in population data, fill them into the Call summary
%along with the corresponding variance frequency in the population data 
if ~isempty(there_idx)
    for j=1:length(there_idx)
        index=there_idx(j);
        Summary(i).SNPs(j).pos=snpsPop(index).pos;
        Summary(i).SNPs(j).varfreq=snpsPop(index).varfreq;
    end
else 
    Summary(i).SNPs=[];
end    


clear j 

%fill mean variance frequency of each sample in population sequencing data
%in CallSummary
if  ~isempty(Summary(i).SNPs)
    Summary(i).MeanFreq=mean([Summary(i).SNPs.varfreq]);
else
    Summary(i).MeanFreq=0;
end
      
%count gaps between found SNPs    
if Summary(i).RatioFound~=1
L=0;
l=1;%running variable
for k=1:length(there)
        if there(k)==0 %if a SNP is missing end 1 to gap length
           L=L+1;
        elseif there(k)~=0 && L~=0 %if a gap ends fill in Summary and reset gap length
            Summary(i).Gaps(l).gaps=L;
            L=0;
            l=l+1;
        else
            continue
        end   
        
        if k==length(there) && there(k)==0 % special case when SNPs are missing at the end of a cluster
            Summary(i).Gaps(l).gaps=L;
        end    
end
clear   k 


%compute mean and std of gapsize
Summary(i).GapMeanPM=[mean([Summary(i).Gaps(:).gaps]), std([Summary(i).Gaps(:).gaps])];
Summary(i).FalseNegatives=snpsMono(~where);
else
    continue
end    


end

%%false positives

%all SNPs occurrent in monoclonal data
allMonos=[MonoData(1).FilterSummary.filtered.pos, MonoData(2).FilterSummary.filtered.pos, ...
    MonoData(3).FilterSummary.filtered.pos, MonoData(4).FilterSummary.filtered.pos];

%find SNPs that are in population SNP/BCFSummary but not in the monoclonal data
[where_fp, there_fp]=ismember([snpsPop(:).pos], allMonos);
fp_pos=~where_fp;


%fill info about false positives in Summary
FalsePositives=snpsPop(fp_pos);
%fill Summary with number of falls negatives and correct calls
NoFalsePositives=size(FalsePositives);
NoFalsePositives=NoFalsePositives(1)

Gaps=[Summary(:).Gaps];
NoFalseNegatives=sum([Gaps.gaps]);
NoCorrectCalls=size([Summary.SNPs]);


%fill CallSummarty with all remaining information 
CallSummary(name_idx).SampleName=SampleData.Sample;
CallSummary(name_idx).NoFalsePositives=NoFalsePositives;
CallSummary(name_idx).FalsePositives=FalsePositives;
CallSummary(name_idx).NoFalseNegatives=NoFalseNegatives;
CallSummary(name_idx).NoCorrectCalls=NoCorrectCalls(2);
CallSummary(name_idx).Info=Summary;


clear Summary NoFalseNegatives NoFalsePositives 
end

%save CallSummary if saveSum=true  
if saveSum
    save(SNPPath + datestr(now, 'yyyymmdd') + "_Mix11_CallSummary.mat", 'CallSummary');
end
