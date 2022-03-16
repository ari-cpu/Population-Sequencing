%% possible figures to visualize the performance of the calling pipelines on SNP level 
%input: CNPSummary, masterlist, changelists

clearvars

path="/home/isabel/shared_Master/Matlab/"; %path to the CNPSummaries 

masterlist=path+"mlBval2Bs166NCe_v1.txt";

changelist=path+"mockfastas/changelistSummary.txt";

SumNames=["20220315_Mix10_CNPSummary.mat" "20220315_Mix11_CNPSummary.mat" "20220315_Mix12_CNPSummary.mat" "20220315_Mix13_CNPSummary.mat"];


%freqs=[summary.CNPSummary.Info(:).varfreq];

%vfhist=histogram(freqs);

%open masterlist
fid=fopen(masterlist);
imp1=textscan(fid, '%f %s %s');
clear fid 
mlpos=imp1{:, 1};
mlref=imp1{2};
mlalt=imp1{3};


%open changelist 
fid=fopen(changelist);
imp2=textscan(fid, '%f %f %f %f %f %f %*f %*f %*f' , 'HeaderLines', 1);
clear fid 

%this list is used to assign each cluster to a known sample
cluster_proben=["Sample 4", "Sample 3", "Sample 3", "Sample 3", "Sample 3", "Sample 2", "Sample 2", "Sample 2"...
    "Sample 1", "Sample 1", "Sample 1"];

%start and end positions from clusters by changelist 
StartCh=imp2{:, 1};
EndCh=imp2{:, 2};

clLen=length(StartCh);
ArtClu(clLen)=struct();

for s=1:length(SumNames) %loop over all CNPSummaries

summary=load(path+SumNames(s));
FoundClu(clLen)=struct();

%fill structure with all positions of SNPs expected to be found in cluster
%from changelist
%also find false positives which lay within the boundary of that cluster
for i=1:clLen
   maskCh=[mlpos(:)]>=StartCh(i) & [mlpos(:)]<=EndCh(i);
   ArtClu(i).clusters=struct('pos', num2cell(mlpos(maskCh)), 'ref',mlref(maskCh), 'alt', mlalt(maskCh));
   ArtClu(i).sample=cluster_proben(i);
   
   FoundClu(i).sample=cluster_proben(i);
   maskfp=[summary.CNPSummary.FPInfo(:).pos]>=StartCh(i) & [summary.CNPSummary.FPInfo(:).pos]<=EndCh(i);
   FoundClu(i).FP=struct('pos', num2cell(summary.CNPSummary.FPInfo(:).pos(maskfp)), 'varfreq', num2cell(summary.CNPSummary.FPInfo(:).vf(maskfp)));
end   

%fill summary with found and missing SNPs within expected clusters
for i=1:clLen
    idx=ismember([summary.CNPSummary.Info(:).pos], [ArtClu(i).clusters(:).pos]);
    idx2=ismember( [ArtClu(i).clusters(:).pos], [summary.CNPSummary.Info(:).pos]);
    
    FoundClu(i).clusters=struct('pos', num2cell(summary.CNPSummary.Info(:).pos(idx)), 'alt', ...
        num2cell(summary.CNPSummary.Info(:).alt(idx)), 'varfreq', num2cell(summary.CNPSummary.Info(:).varfreq(idx)));
    
    FoundClu(i).missing=struct('pos', num2cell([ArtClu(i).clusters(~idx2).pos]));
    
    %trying to compute positions where cluster is split due to 5
    %consectucive missing SNPs but it does not work properly yet
    %aaaaaaaaaaaaaah
    spidx=~idx2;
    splitpos=find((spidx.*circshift(spidx, 1).*circshift(spidx, 2).*circshift(spidx, -1).*circshift(spidx, -2)));
    splitposrot=splitpos-circshift(splitpos, 1);
   
    if ~isempty(find(circshift(splitposrot, -1)>1))
    
    spMask(1, :)=find(circshift(splitposrot, -1)>1);
    spMask(2, :)=find(circshift(splitposrot, -1)>1)+1;
    sppos=splitpos(spMask);
    FoundClu(i).SplitPositions=[FoundClu(i).clusters(sppos).pos]
   
    elseif isempty(splitposrot)
       
       FoundClu(i).SplitPositions=[];
    %else
       
     %  sppos=splitpos(1); 
      % if FoundClu(i).clusters(sppos).pos >= min([FoundClu(i).clusters(:).pos])
      % FoundClu(i).SplitPositions=[FoundClu(i).clusters(sppos).pos];
      % else     
      % FoundClu(i).SplitPositions=[];
      % end
           
       
    end
    
    %fill the percebtage of SNPs found in each cluster in structure
    FoundClu(i).ratioFound=length([FoundClu(i).clusters(:).pos])/length([ArtClu(i).clusters(:).pos]);
   
    
    
end

    FoundSummary(s).MixData=FoundClu;
    FoundSummary(s).MixName=SumNames(s);
    clear FoundClu
    
end
    

mix=1; %index of mix you want to visualize
FoundClu=FoundSummary(mix).MixData;
%%  plotssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

%frequency histograms for each cluster seperately


for id=1:clLen
subplot(3, 4, id);
vf1=histogram([FoundClu(id).clusters(:).varfreq]);
title(FoundClu(id).sample);
xlabel("variance frequency");
if ~isempty([FoundClu(id).clusters(:).varfreq])
meanvf=mean([FoundClu(id).clusters(:).varfreq]);
stdvf=std([FoundClu(id).clusters(:).varfreq]);
legend(string(meanvf) + newline + string(stdvf), 'Location','northwest');
end
end


%posart=[1, 2, 3, 7, 8, 9, 13, 14, 15, 19, 20];
%posfound=[4, 5, 6, 10, 11, 12, 16, 17, 18, 22, 23];

%% visualization of clusters
%plot depicting each position in each cluster color coded: green -> found
%SNP, red -> missing SNP, blue -> false positive
%very slow for some reason
for id=1:clLen
 
    subplot(4, 3, id);
    foundclu=[FoundClu(id).clusters(:).pos];
    foundcluy=zeros(length(foundclu));
    missclu=[FoundClu(id).missing(:).pos];
    misscluy=zeros(length(missclu));
    if ~isempty(FoundClu(id).SplitPositions)
        for i=1:length(FoundClu(id).SplitPositions')
        hold on;
        xline(FoundClu(id).SplitPositions(i))
        end
    end    
    fpclu=[FoundClu(id).FP(:).pos];
    fpcluy=zeros(length(fpclu));
    plot(foundclu, foundcluy, 'g|');
    hold on;
    plot(missclu, misscluy, 'r|');
    hold on;
    plot(fpclu, fpcluy, 'b|');
    g=gca;
    ylim([0, 0.01]);
    title(FoundClu(id).sample);
    %xlim=get(g, 'Xlim');
    %xlabel('position on genome', 'Position', [xlim(1) -0.1]);
    set(g, 'box', 'off');
    g.YAxis.Visible='off';
    g.XRuler.TickLabelGapOffset=10;
    g.XAxis.TickLength=[0 0];
    
    clear foundclu foundcluy missclu misscluy
end    

%% visualization with colormap 
% plot depicting each position in the clusters including found SNPs and
% false positives
%the colormap corresponds to the found variant frequency at that position 
%still lacks legend etc
samples=["Sample 1", "Sample 2", "Sample 3", "Sample 4"];
colormap(cool);
tiledlayout('flow')
for id=4:4%length(samples)
    
    sampleClu=FoundClu(find([FoundClu(:).sample]==samples(id)));
    len=length(sampleClu);
    vfs=[];
    
    for j=1:len
        vfs=[vfs [sampleClu(j).clusters(:).varfreq] [sampleClu(j).FP(:).varfreq]];
    end    
    
    
    for a=1:len
    
  
    foundclu=[sampleClu(a).clusters(:).pos];
    foundcluy=zeros(1, length(foundclu));
    fpclu=[sampleClu(a).FP(:).pos];
    fpcluy=zeros(1, length(fpclu));
    
    c1=[sampleClu(a).clusters(:).varfreq];
    c2=[sampleClu(a).FP(:).varfreq];
    nexttile
    scatter(foundclu, foundcluy,100,c1, '|');
    hold on;
    scatter(fpclu, fpcluy,100,c2, 'x');
    colorbar
    caxis([min(vfs) max(vfs)])
    title(samples(id));
    ylim([0 0.01]);
    p=gca;
    set(p, 'box', 'off');
    p.YAxis.Visible='off';
    p.XAxis.TickLength=[0 0];
    p.XRuler.TickLabelGapOffset=8;
    end
end




%% old plot, not really useful atm but I may want to keep it for later

colormap(winter);
  [C, I]=sort([FoundClu(a).clusters(:).varfreq]);
        ei=max(C)-min(C);
        Cscale=[];


        for i=1:length(C)
            Cscale=[Cscale (max(C)-C(i))*(1/ei)];
    
        end

    Cscale=1-Cscale;
    zmap=linspace(min(Cscale), max(Cscale), length(cmap));
    tickmap=linspace(min(C), max(C), 6);
for id=1:clLen
    subplot(4, 3, id)
    i=I(id);
    foundclu=[FoundClu(i).clusters(:).pos];
    foundcluy=zeros(length(foundclu));
    color=interp1(zmap, cmap, Cscale(id));
   
    plot(foundclu, foundcluy, '|', 'color', color);
    cb=colorbar();
    cb.Ticks=linspace(0, 1, 6);
    cb.TickLabels=(num2cell(tickmap));
    title(FoundClu(i).sample);
    p=gca;
    ylim([0 0.01]);
    set(p, 'box', 'off');
    p.YAxis.Visible='off';
    p.XAxis.TickLength=[0 0];
    legend('detected ratio:'+string(C(id)));
   
    
   
end


%% dor plot containing info about all four mixes in comparison 
%you can vary betwen yAxis = mean frequency found for cluster ([data(:).meanfr])
%or YAxis= ratio of SNPs correctly found in each cluster ([data(:).ratioFound])

sampleName=["Sample 1" "Sample 2" "Sample 3" "Sample 4"];
Samplefreq=[ 0.533 0.267 0.133 0.067];
color=[1 2 3 4];
colormap(jet)
for u=1:4
  data=FoundSummary(u).MixData;
  sampleFreqShift=circshift([Samplefreq], u-1);
  for t=1:length(data)
      idx3=find(sampleName==[data(t).sample]);
      data(t).sampleFreq=sampleFreqShift(idx3);
      data(t).c=color(idx3);
      data(t).meanfr=mean([data(t).clusters(:).varfreq]);
      data(t).NoSNPs=length(data(t).clusters)
  end
  scatter([data(:).sampleFreq], [data(:).ratioFound],100,[data(:).c], '*');
  colorbar;
  hold on;
end    
