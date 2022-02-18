clearvars

                     %                      %
                    % %    Define input    % %
                     %                      %

donor = "Bval" 
sample = "Bsub_MockReads_Mix3_changeAll" 
saveSum=true;


SNPPath = "/home/isabel/Documents/shared_Master/Matlab/";
artefacts=SNPPath+"Artefacts_Bs1662BsNCe.vcf";

SNPName{1} = "20220217_SNPSummary.mat";
%SNPName{2} = "20211124_comp_SNPSummary.mat";
%SNPName{2} = "20220210_comp_SNPSummary.mat"; 
SNPName{2} = "20220217_comp_SNPSummary.mat"; 


PopData=load(SNPPath+SNPName{1});
MonoData=load(SNPPath+SNPName{2});

PopData=PopData.SNPSummary;
MonoData=MonoData.SNPSummary;

name_idx=[PopData.Sample]==sample;

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

%% comparison with the Monoclonal Data 
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
%check if positions of monosnps are found on pool data
[where, there]=ismember([snpsMono(:).pos], [snpsPop(:).pos]);
there_idx=there(where);

%if isempty(there_idx)
%  continue
%end  

%check for faulty snps
checkSNPs=[snpsMono(where).alt]==[snpsPop(there_idx).alt];

if ~isempty(nonzeros(checkSNPs))
    Summary(i).NoFalseSNPs=0;
    Summary(i).FalseSNPS=[];
else
    Summary(i).NoFalseSNPs=length(nonzeros(checkSNPs(~(nonzeros(checkSNPs)))));
    Summary(i).FalseSNPs=snpsPop(~checkSNPs);
end

%compute ratio of snps that is found in pool data
Summary(i).RatioFound=length(nonzeros(where))/length([snpsMono(:).pos]);
%Summary(i).RatioInTheSample=length(nonzeros(where))/length(snpsPopPos);


  

%list all found SNPs in pool data
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

if  ~isempty(Summary(i).SNPs)
    Summary(i).MeanFreq=mean([Summary(i).SNPs.varfreq]);
else
    Summary(i).MeanFreq=0;
end
      
%count gaps between found SNPs    
if Summary(i).RatioFound~=1
L=0;
l=1;
for k=1:length(there)
        if there(k)==0 
           L=L+1;
        elseif there(k)~=0 && L~=0
            Summary(i).Gaps(l).gaps=L;
            L=0;
            l=l+1;
        else
            continue
        end   
        
        if k==length(there) && there(k)==0
            Summary(i).Gaps(l).gaps=L;
        end    
end
clear   k 



%mean and std of gapsize
Summary(i).GapMeanPM=[mean([Summary(i).Gaps(:).gaps]), std([Summary(i).Gaps(:).gaps])];
Summary(i).FalseNegatives=snpsMono(~where);
else
    continue
end    


end

%%false positives

allMonos=[MonoData(1).FilterSummary.filtered.pos, MonoData(2).FilterSummary.filtered.pos, ...
    MonoData(3).FilterSummary.filtered.pos, MonoData(4).FilterSummary.filtered.pos];

[where_fp, there_fp]=ismember([snpsPop(:).pos], allMonos);
fp_pos=~where_fp;


FalsePositives=snpsPop(fp_pos);
NoFalsePositives=size(FalsePositives);
NoFalsePositives=NoFalsePositives(1)

Gaps=[Summary(:).Gaps];
NoFalseNegatives=sum([Gaps.gaps]);
NoCorrectCalls=size([Summary.SNPs]);

CallSummary.SampleName=sample;
CallSummary.NoFalsePositives=NoFalsePositives;
CallSummary.FalsePositives=FalsePositives;
CallSummary.NoFalseNegatives=NoFalseNegatives;
CallSummary.NoCorrectCalls=NoCorrectCalls(2);
CallSummary.Info=Summary;

if saveSum
    save(SNPPath + datestr(now, 'yyyymmdd') + "_CallSummary.mat", 'CallSummary');
end    

