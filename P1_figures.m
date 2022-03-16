%% just a bunch of possible figures to visualize the performance of the calling pipelines on SNP level 
%input: CallSummary
clearvars

donor="Bval"
sample="Bsub_MockGenome_Mix_changeAll_bcfmp"

SumPath = "/home/isabel/shared_Master/Matlab/";


CallNameBCFtools="20220314_Mix10_BCFcall_CallSummary.mat" %Summary generated ftom BCFtools call function
CallNameBCF="20220314_Mix_DP_CallSummary_ml.mat" %Summary generated from custom filtering of bcf.vcf files
CallNameGATK="20220314_Mix_gatk_CallSummary.mat" % Summary generated from GATK Haplotype Caller

CallBCF=load(SumPath+CallNameBCF);
CallGATK=load(SumPath+CallNameBCFtools);

InfoBCF=CallBCF.CallSummary.Info;
InfoGATK=CallGATK.CallSummary.Info;


ratio=[0.5, 0.25, 0.125, 0.0625];


%% generate arrays with the number of false nagtives and false positives counting into each ratio categrory
fnArrBCF=[zeros(length(InfoBCF(1).FalseNegatives), 1)+ratio(1); zeros(length(InfoBCF(2).FalseNegatives), 1)+ratio(2);...
   zeros(length(InfoBCF(3).FalseNegatives), 1)+ratio(3); zeros(length(InfoBCF(4).FalseNegatives), 1)+ratio(4)];
fnArrGATK=[zeros(length(InfoGATK(1).FalseNegatives), 1)+ratio(1); zeros(length(InfoGATK(2).FalseNegatives), 1)+ratio(2);...
   zeros(length(InfoGATK(3).FalseNegatives), 1)+ratio(3); zeros(length(InfoGATK(4).FalseNegatives), 1)+ratio(4)];

fpBCF=[CallBCF.CallSummary(:).FalsePositives.varfreq];
fpGATK=[CallGATK.CallSummary(:).FalsePositives.varfreq];

%for the false positives, the ratio they belong to is just with in half of
%the gaps between ratios
%no false positive is discarded
fpcountBCF=flip([length(fpBCF(fpBCF<=0.09375)), length(fpBCF(fpBCF>0.09375 & fpBCF<=0.1875)),...
    length(fpBCF(fpBCF>0.1875 & fpBCF<=0.375)), length(fpBCF(fpBCF>0.375))]);

fpcountGATK=flip([length(fpGATK(fpGATK<=0.09375)), length(fpGATK(fpGATK>0.09375 & fpGATK<=0.1875)),...
    length(fpGATK(fpGATK>0.1875 & fpGATK<=0.375)), length(fpGATK(fpGATK>0.375))]);

fpratioBCF=fpcountBCF./[InfoBCF(:).numberOfSNPs];
fpratioGATK=fpcountGATK./[InfoGATK(:).numberOfSNPs];
%% bar plot comparing the paerformance of two callers by the number of false negatives, false positives and correct calls 
tiledlayout(1, 1);
ax=nexttile;
foundratiosBCF=[[InfoBCF(:).RatioFound]; 1-[InfoBCF(:).RatioFound]; fpratioBCF];
foundratiosGATK=[[InfoGATK(:).RatioFound]; 1-[InfoGATK(:).RatioFound]; fpratioGATK];
stack1=bar([0.8, 1.8, 2.8, 3.8], foundratiosBCF, 0.3,'stacked');
stack1(1).FaceColor='green';
stack1(2).FaceColor='magenta';
stack1(3).FaceColor='cyan';
hold on;
stack2=bar([1.2, 2.2, 3.2, 4.2], foundratiosGATK,0.3, 'stacked');
stack2(1).FaceColor='green';
stack2(2).FaceColor='magenta';
stack2(3).FaceColor='cyan';

ticks={'custom', '0.5', 'BCFtools', 'custom', '0.25', 'BCFtools', 'custom', '0.125', 'BCFtools', 'custom', '0.0625', 'BCFtools'};
xticks([0.8 1 1.2 1.8 2 2.2 2.8 3 3.2 3.8 4 4.2]);
xticklabels(ticks);
ylim([0, 1.5]);
legend('{SNP found}', '{SNP not found}', '{False Positive}');
xlabel('variant frequency');

%% histogram comparing the false positves varinat frequency of two callers
edges=[0:0.05:0.4];
hist1=histogram(fpBCF, edges);
hold on;
hist2=histogram(fpGATK, edges);
%xlim([0, 0.1]);
%set(gca, 'Yscale', 'log');
legend('{BCF.vcf cutoff:20}', '{GATK}');
xlabel('variant frequency');
title('variant frequencies of false positives');


%% bcftools call performance 
%line plot with the fraction mixed into the mock genome vs the fraction of
%SNPs found by BCFtools call 
fraction=[0.1, 0.15, 0.2, 0.25, 0.5, 1];
found1=[0, 0, 0.17, 0.61, 0.99, 1];
found2=[0, 0, 0.07, 0.56, 0.98, 1];
found3=[0, 0, 0.07, 0.54, 0.98, 1];
found4=[0, 0, 0.07, 0.59, 0.99, 1];
foundMean=[0, 0, (found1(3)+found2(3)+found3(3)+found4(3))/4, (found1(4)+found2(4)+found3(4)+found4(4))/4, ...
    (found1(5)+found2(5)+found3(5)+found4(5))/4, 1];
%found1=[1, 1, 1, 1, 1, 1];
%found2=[0.98, 1, 1, 1, 1, 1];
%found3=[0.98, 0.98, 0.98, 0.98, 0.98, 0.98];
%found4=[1, 1, 1, 1, 1, 1];

p1=plot(fraction, found1, 'mo','DisplayName', 'Mockgenome 1');
hold on;
p2=plot(fraction, found2, 'c*','DisplayName', 'Mockgenome 2');
hold on;
p3=plot(fraction, found3, 'gd','DisplayName', 'Mockgenome 3');
hold on;
p4=plot(fraction, found4, 'b+','DisplayName', 'Mockgenome 4');
hold on; 
p5=plot(fraction, foundMean, 'k--', 'DisplayName', 'Mean value');
h=gca;
set(h, 'box', 'off');
set(h, 'TickLength', [0.005 0.005]);

legend('Location', 'northwest');
xlabel('fraction in the sample');
ylabel('fraction found by bcftools call');
