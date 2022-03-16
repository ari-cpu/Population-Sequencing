%comparing the SNPs found  by monoclonal ccalling to SNPs expected from the
%masterlist and changelists
%input: SNPSummary from monoclonal data, changelists, masterlist

clearvars

                     %                      %
                    % %    Define input    % %
                     %                      %

donor = "Bval" 
sampleName="Bsub_MockReads_change"


SNPPath = "/home/isabel/shared_Master/Matlab/"; 
artefacts=SNPPath+"Artefacts_Bs1662BsNCe.vcf";
masterlist=SNPPath+"mlBval2Bs166NCe_v1.txt";
change=SNPPath+"mockfastas/changelist";
%fasta=SNPPath+"mock_reads/"+SampleName;
%SNPName="20211124_comp_SNPSummary.mat";
SNPName="20220314_newMixes_comp_SNPSummary.mat" %name of SNPSummary

MonoData=load(SNPPath+SNPName); %load data from SNPSummary

% load list of artefacts
if ~isempty(artefacts)
fid = fopen(artefacts);
imp = textscan(fid, '%s %f %s %s %s %f %s %s %s %s', 'HeaderLines', 33);
artefact.pos = imp{2};
artefact.alt = imp{5};
fclose(fid);
clear imp fid;
end 

%open Masterlist

fid=fopen(masterlist);
imp=textscan(fid, '%f %s %s');
%ml=cell2struct(imp, {'pos', 'ref', 'alt'}, 2);
ml=struct('pos', num2cell(imp{1, 1}), 'ref', imp{1, 2}, 'alt', imp{1, 3});
fclose(fid);
clear imp fid;






G=4;

%initialize structure for info on performance of calling
CompSummary(G)=struct();
%%
% loop over all mockgenomes that are to be tested
for i=1:G
    changelist=change+int2str(i)+".txt";
    Mono=MonoData.SNPSummary(i).FilterSummary.filtered; %%%chnagethisback!!!!! 
    
    %filter out artefacts if there are any
    for a = 1 : length(artefact.pos)
    idx = find([Mono(:).pos] == artefact.pos(a));
    if ~isempty(idx) %&& strcmp(artefact.alt(a), SampleData.FilterSummary.filtered(idx).alt)
       Mono(idx)=[];
    end   
    end 
    clear idx a
  
    %open changelist of the mockgenome
    fid=fopen(changelist);
    imp=textscan(fid, '%f %f %f %f %*f %*f %*f %*f %*f %*f', 'Headerlines', 1);
    cl=cell2struct(imp, {'StartSbj', 'EndSbj', 'StartQry', 'EndQry'}, 2);
    fclose(fid);
    clear imp fid;
    
    CompSummary(i).MockName="changelist"+int2str(i);
    CompSummary(i).NoSNPs=length(Mono);
    CompSummary(i).Info=struct();
    
   clLen=length(cl.StartSbj);
   maskml=zeros(length([ml(:).pos]), clLen);
   
  
   % create seperate mask for each piece of genome in this changelist for masterlist
   for p=1:clLen
       maskml(:, p)= [ml(:).pos]<cl.StartSbj(p) | [ml(:).pos]>cl.EndSbj(p); 
   end    
   
   %Combine masks for all positions in that changelist
   maskml=~logical(prod(maskml, 2));
   
   %list of all SNPs expected from the changelists 
   ml_i=ml(maskml);
   
   %% missing SNPs, faulty SNPs, false positives  
  
   %find missing SNPs
   mlthere=ismember([ml_i(:).pos], [Mono(:).pos]);
   CompSummary(i).Info.MissingSNPs=ml_i(~mlthere);
   CompSummary(i).NoMissingSNPs=length(nonzeros(~mlthere));
   
   %false positives within the clusters
   [monothere, monoloc]=ismember([Mono(:).pos], [ml_i(:).pos]);
   monoloc=monoloc(monothere)';
   monothere=monothere';
   CompSummary(i).Info.FalsePositives=Mono(~monothere);
   CompSummary(i).NoFalsePositives=length(nonzeros(~monothere));
   
   % faulty SNPs
 
   Mono=Mono(monothere);
   ml_i=ml_i(monoloc);
   
   faults=find([Mono(:).alt]~=[ml_i(:).alt]);
   CompSummary(i).Info.FaultySNPs=struct('pos', num2cell([Mono(faults).pos]), 'ref', num2cell([Mono(faults).ref]), 'altMock', ...
   num2cell([Mono(faults).alt]), 'altMl', num2cell([ml_i(faults).alt]));
   CompSummary(i).NoFaultySNPs=length(faults);
   
   
 

   
    
end

    
