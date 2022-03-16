%input: SNPSummary or BCFSummary, masterlist, changelists
%generates a CallSummary containing found SNPs compared to expected SNPs,
%false positives etc
%basis is the comparison with expected SNPs from the masterlist

clearvars

                     %                      %
                    % %    Define input    % %
                     %                      %

donor = "Bval" 
saveSum=true; %save the CallSummary ? 


SNPPath = "/home/isabel/shared_Master/Matlab/"; %path where SNP/BCFSummary is found
chpath=  "/home/isabel/shared_Master/Matlab/mockfastas/"; %path where changelists are found
artefacts=SNPPath+"Artefacts_Bs1662BsNCe.vcf";
masterlist=SNPPath+"mlBval2Bs166NCe_v1.txt";
SNPName{1} = "20220315_Mix13_BCFSummary.mat"; %name of the Summary you want to analyze
changelist=chpath+"changelist";
%SNPName{1}="20220310_Mix1_BCFSummary.mat";


PopData=load(SNPPath+SNPName{1}); %load Summary of population sequencing SNPs
PopData=PopData.SNPSummary;


fid=fopen(masterlist); %open masterlist
imp=textscan(fid, '%f %s %s');
%ml=cell2struct(imp, {'pos', 'ref', 'alt'}, 2);
ml=struct('pos', num2cell(imp{1, 1}), 'ref', imp{1, 2}, 'alt', imp{1, 3});
fclose(fid);
clear imp fid;



mockChanges(4)=struct();
cluster_proben=["Sample 1", "Sample 2", "Sample 3", "Sample 4"];
cluster_freq=circshift([0.533 0.267 0.133 0.067], 3); % these are the rations of each sample in the artificial mix 


for g=1:4 % create a structure containing the start and end positions of the artificial clusters
    fid=fopen(changelist+string(g)+".txt");
    imp2=textscan(fid, '%f %f %f %f %f %f %*f %*f %*f' , 'HeaderLines', 1);
    clear fid 
    
    mockChanges(g).sample=cluster_proben(g);
    
    mockChanges(g).StartCh=imp2{:, 1};
    mockChanges(g).EndCh=imp2{:, 2};

 
end



ArtClu(4)=struct();

for i=1:4 %fill a structure with the expected SNPs in the artificial clusters
    StartCh=mockChanges(i).StartCh;
    EndCh= mockChanges(i).EndCh;
    maskCh=[];
   
    for n=1:length(StartCh)
    maskCh=[maskCh; [ml(:).pos]>=StartCh(n) & [ml(:).pos]<=EndCh(n)];
    end
    mlMask=logical(sum(maskCh, 1));
    
    
    ArtClu(i).clusters=struct('pos', num2cell([ml(mlMask).pos]), 'ref',num2cell([ml(mlMask).ref]), 'alt', num2cell([ml(mlMask).alt]));
 
    ArtClu(i).sample=cluster_proben(i);
    ArtClu(i).sampleFreq=cluster_freq(i);
    clear n mlMask maskch
   
end   




%name_idx=[PopData.Sample]==sample;
i=size(PopData);

for name_idx=1:i(2) % loop over all SNP/BCFSummaries

SampleData=PopData(name_idx);

snpsPop=[SampleData.FilterSummary.filtered]; % filter summary of current sample


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




%% comparison with the SNPs expected from masterlist
for i=1:length(ArtClu) %loop over all samples mixed into the populartion

snpsArt=ArtClu(i).clusters;




Summary(i).numberOfSNPs=length([snpsArt(:).pos]); %expected number of SNPs

[where, there]=ismember([snpsArt(:).pos], [snpsPop(:).pos]); %check which of the expected SNPs are found in the population data
there_idx=there(where);


%check for faulty snps (right position bur wrong allele)
checkSNPs=[snpsArt(where).alt]==[snpsPop(there_idx).alt];

if  isempty(nonzeros(checkSNPs))
    Summary(i).NoFalseSNPs=0;
    Summary(i).FalseSNPS=[];
else
    Summary(i).NoFalseSNPs=length(nonzeros(checkSNPs(~(nonzeros(checkSNPs)))));
    Summary(i).FalseSNPs=snpsPop(~checkSNPs);
end

%compute ratio of expected snps that is found in population data
Summary(i).RatioFound=length(nonzeros(where))/length([snpsArt(:).pos]);



  

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

%fill expected variance frequency for each sample in CallSummary
Summary(i).sampleFreq=ArtClu(i).sampleFreq;

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
l=1; %running variable
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
Summary(i).FalseNegatives=snpsArt(~where);
else
    continue
end    


end

%%false positives

%all SNPs occurrent in all clusters
allArt=[ArtClu(1).clusters.pos, ArtClu(2).clusters.pos, ArtClu(3).clusters.pos, ArtClu(4).clusters.pos];

%find SNPs that are in SNP/BCFSummary but not in any expected cluster
[where_fp, there_fp]=ismember([snpsPop(:).pos], allArt);
fp_pos=~where_fp;

%fill info about false positives in Summary
FalsePositives=snpsPop(fp_pos);
NoFalsePositives=size(FalsePositives);
NoFalsePositives=NoFalsePositives(1)

Gaps=[Summary(:).Gaps];
%fill Summary with number of falls negatives and correct calls
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
    save(SNPPath + datestr(now, 'yyyymmdd') + "_Mix13_CallSummary_ml.mat", 'CallSummary');
end