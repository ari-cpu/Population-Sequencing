% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                     VariantFiltering.m                              % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clearvars

             % %                                  % %
            % % %   Specify your samples ...     % % %
             % %                                  % %

filepath = "/home/isabel/shared_Master/Matlab/bcfcall/mock/";    % Where are the unfiltered variant lists 
                                %(*_bcfcall.vcf) and where do you 
                                % want to save the SNPSummary?
%filepath = path+"pipeline_output_"+mode+"\"; % + "RAW\";
path ="/home/isabel/shared_Master/Matlab/"

filesuffix = "_bcfcall.vcf";

savepath = path; % + "Filtered\";
% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

sampleType = "customized";    % Using the standard sample naming (under our 
                            % convention), choose sampleType = "standard";
                            % not standard naming but cycle included in the
                            % name can use "standard" as well;
                            % To filter other files (cycle not in name), 
                            % set sampleType = "customized"; Note, that for 
                            % these files LAfiltering is not possible !

% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %
% % % --- --- --- --- ---  sampleType == "standard" - --- --- --- --- --- % % %
% % % --- --- In this sampleType, the naming of a sample follows: --- --- % % %
% % % Prefix(evolExp) + sampleNo + corrCycle + sampleSuffix + filesuffix  % % %
% % % e.g. Wns        +    01    +    10     +   _2Donor    + _bcfcall.vcf% % %
% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

if strcmp(sampleType, "standard")
    % % evolExp == samplePrefix
    evolExp = "Vns";      

    % Define sample number sampleNo{x} and the corresponding cycles corrCycles{x}:
    sampleNo{1} = [ "02", "03", "07", "08"];
    corrCycles{1} = [10];

    sampleNo{2} = ["01", "04", "05", "06", "09", "10", "11", "12", "13", "14" ];
    corrCycles{2} = [10 20];

    sampleSuffix = "";  % You can specify a sampleSuffix if necessary, default:
                        % sampleSuffix = "";
    %sampleSuffix = "_2Donor";
    
% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %
% % % --- --- --- --- --- sampleType == "customized"  --- --- --- --- --- % % %
% % % In the "customized" - sampleType, the entire sample names need to - % % %
% % % be fill in the sampleNames string array, leave the rest untouched - % % % 
% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

elseif strcmp(sampleType, "customized")
    
    %sampleNames = ["LibSCBval1_2NCe" "LibSCBval7_2NCe" "LibSCBval9_2NCe" "LibSCBval17_2NCe" "LibSCBval89_2NCe"];
    %sampleNames = ["LibBval_pop_equal_1" "LibBval_pop_equal_2" "LibBval_pop_equal_3" "LibBval_pop_B5_H5" "LibBval_pop_A9_B5" "LibBval_pop_H5_A1"];
    sampleNames = ["Bsub_MockReads_newMixes_change1_5" "Bsub_MockReads_newMixes_change2_5" "Bsub_MockReads_newMixes_change3_5" "Bsub_MockReads_newMixes_change4_5"];
    %sampleNames=["Bsub_MockReads_Mix10_changeAll_bcfmp"]
    
    corrCycles{1} = [0];
    sampleNo{1} = sampleNames; evolExp = "";
    
else
    error("Please, clarify you sampleType!");
end 

% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %
% % % --- ---  What filtering parameters would you like to set ?  --- --- % % %
% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

minTotalCount = 50;     % What is your total count cut off ? (default: 50)
totalCountsMode = "DP"; % Do you like the raw read depth DP or the sum over 
                        % all allelic depths AD be taken as your total 
                        % counts ? ("DP" / "AD")
minQual = 50;        % What is the minimum phred quality score? (default: 50)
unambiFreq = 90;     % [in %]; When is a variant unambigious ? (default : 90)

% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

% Do you want to save the SNPSummary to the savepath ?
saveSNPSummary = true;

% % % --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- % % %

  %                                                                  %
 % %                                                                % %
% % %   Now, just lean back and let me do the rest for you ...     % % %
 % %                                                                % %
  %                                                                  %
  
SNPSummary = struct('IndvMutList',[], 'Sample',[], 'LAfiltered', [], 'FilterSummary', []);

% allSamples is a list of all samples (numbers), the corresponding cycles 
% can be found in allCycles:
allSamples = [sampleNo{:}];
allCycles = [];
for i = 1 : length(sampleNo)
   allCycles = [allCycles; repmat(corrCycles(i), length(sampleNo{i}),1)];
end

% Loop over all samples
allCount = 1;
for sampleCount = 1 : length(allSamples)
    
    % Pick current sample
    sampleCurr = evolExp + allSamples(sampleCount);
    % What cycles do we provide for this sample ?
    cycles = allCycles{sampleCount};

    if strcmp(sampleType, "standard")
        fileName = cellstr(strcat(filepath, sampleCurr, string(num2str(cycles','%02.f')), sampleSuffix, filesuffix))';
    elseif strcmp(sampleType, "customized")
        fileName = cellstr(strcat(filepath, sampleCurr, filesuffix))';
    end
    
    FilterSummary = struct(...
            "cycle",num2cell(cycles), "filename", fileName,  ...
            "filterDate", datestr(now, 'dd/mm/yyyy'), "totalCountMode", totalCountsMode, ...
            "minTotalCounts", minTotalCount, "minQual", minQual,...
            "unambigiousFreq", unambiFreq, "LAfiltered", [], ...
            "unfiltered",cell(size(cycles)), ...
            "filtered",cell(size(cycles))...
        );

    % Go through all cycles, starting from the latest one: 
    for cycleBackCounter = 1 : length(cycles)
        
        %% Read data from bcf-call
        fprintf("Read bcf-call '%s'\n", FilterSummary(end-cycleBackCounter+1).filename);

        fid = fopen(FilterSummary(end-cycleBackCounter+1).filename);
        imp = textscan(fid,'%*s %d . %s %s %f . %s %*s %s', 'CommentStyle','#'); % Columns with %* will not be read in
        fclose(fid);

        % In the struct snp, you will now write Position, Reference, Alternate, 
        % Info and Counts :
        snp = cell2struct(horzcat(num2cell(imp{1}),imp{2:3}, num2cell(imp{4}), imp{5:end}), ...
            ["pos","ref","alt","qual","info","counts"],2);
        clear imp fid;

        %% Convert bcf-call data to a useable data structure

        fprintf("\tConvert data...\n");

        splittedCounts = split({snp.counts}',":"); %just splitting last bcfcall-column at :
        altsArr = {snp.alt};    % all alternates
        altCounts = count({snp.alt}',",") + 1;  % count number of alternates

        % Preallocate memory for countArrs and altArrs
        countArrs = cell(size(snp));
        altArrs = cell(size(snp));
        % Splitting splittedCounts and altsArr:
        for alts = 1:3
            maskPos = altCounts == alts;
            countArrs(maskPos) = num2cell(str2double(split(splittedCounts(maskPos,end),",",2)),2);
            altArrs(maskPos) = mat2cell(split(altsArr(maskPos)',",",2),ones(sum(maskPos),1));
        end
        clear splittedCounts altsArr altCounts alts mask;

        totalCounts = cellfun(@(x) sum(x), countArrs);

        % Now, decide for the alt with the higher frequency (more counts)
        
        altAllele = cellfun(@(x) outputs2array(@()max(x(2:end)),[1 2]), countArrs, 'UniformOutput',false);
        altAllele = vertcat(altAllele{:});
        altAlleleChar = cellfun(@(alt,idx) alt{idx}, altArrs, num2cell(altAllele(:,2)), 'UniformOutput', false);
        altFrequency = altAllele(:,1) ./ totalCounts;

        clear altArrs altAllele;
        
       if strcmp(totalCountsMode, 'DP')
            % Extract DP info
            cellY = cellfun(@(s) strsplit(s, 'DP='), {snp.info}', 'UniformOutput', false);
            cellZ = cellfun(@(s) strsplit(s{2}, ';'), cellY, 'UniformOutput', false);
            cellDP = cellfun(@(s) str2double(s{1}), cellZ, 'UniformOutput', false);
            DP = cell2mat(cellDP);
            clear cell*
        else
        end
        
        %% Creating filters ...
        fprintf("\tCreate filter...\n");

        filterIndels = ~contains({snp.info},"INDEL");
        
        if strcmp(totalCountsMode, 'AD')
            filterTotalCounts = totalCounts >= minTotalCount;
        elseif strcmp(totalCountsMode, 'DP')
            filterTotalCounts = DP >= minTotalCount;
        else
            error("Your total count filtering fails! Please set your totalCountsMode to AD or DP.");
        end
        
        filterQual = [snp.qual]' >= minQual; 

        filterUnambi = altFrequency >= unambiFreq / 100; 

        if  cycleBackCounter > 1
            LAfilt = true;
            % 1. Which snps do we find in the following step (checking position)?? --> mask/idx
            [maskPos, idx] = ismember([snp.pos]', [FilterSummary(end-cycleBackCounter+2).filtered.pos]); % mask of unfiltered current snps that are present in the following cycle
            % 2. Check allele:
            prevAlts = horzcat({FilterSummary(end-cycleBackCounter+2).filtered.alt}, {''}); %Add an empty char to variants of following cycle why??
            idx(idx == 0) = length(prevAlts); %Alle Positionen die in idx = 0 waren werden = length(prevAlts) gesetzt
            maskAlt = strcmp(altAlleleChar,prevAlts(idx)');
            % 3. Are they unambigious the next timestep ? 
            maskAmbiCy = vertcat([FilterSummary(end-cycleBackCounter+2).filtered.varfreq]' >= unambiFreq / 100 , 0);% Which alleles are unambigous the following cycle
            maskAmbiCy = maskAmbiCy(idx);
            
            % If a variant is ambigious but(!) unambigious in the
            % following cycle(1.+2.+3.)(-> "Look ahead"), it will be kept:
            filterLookAhead = ~filterUnambi & maskPos & maskAlt & maskAmbiCy;
            
            % Now combining the ambigiousity filters:
            % All snps that are unambigious or ! ambigious but unambigious 
            % in the following time step! will be kept:
            filterComb = filterUnambi | filterLookAhead; 
            
            filterAll = filterIndels' & filterTotalCounts & filterQual & filterComb;

        else  % only for cycleBackCounter = 1 , so the last time step 
            LAfilt = false; % no LookAheadFiltering !
            filterAll = filterIndels' & filterTotalCounts & filterQual & filterUnambi; 
            filterLookAhead = zeros(length(filterAll), 1);            
        end

        %% Write FilterSummary 

        % Write unfiltered data in FilterSummary
        FilterSummary(end - cycleBackCounter + 1).unfiltered = ...
            struct("pos",{snp.pos}', "DP", num2cell(DP), "ref",{snp.ref}', "alt",{snp.alt}', ...
            "varcounts",countArrs,"varfreq", num2cell(altFrequency), "QUAL", {snp.qual}', "filterIndels",num2cell(filterIndels)', "filterQual", num2cell(filterQual), ...
            "filterDP",num2cell(filterTotalCounts), "filterUnambi",num2cell(filterUnambi), ...
            "filterLA", num2cell(filterLookAhead),  "filterAll",num2cell(filterAll)); 

        % Filter data
        snpFilt = snp(filterAll);

        % Write filterd data in FilterSummary
        if ~isempty(snpFilt)
            FilterSummary(end-cycleBackCounter+1).filtered = ...
                struct("pos",{snpFilt.pos}', "ref",{snpFilt.ref}', "alt",altAlleleChar(filterAll), ...
                "varfreq",num2cell(altFrequency(filterAll)));
        else 
            FilterSummary(end-cycleBackCounter+1).filtered = ...
                struct("pos",[], "ref",[], "alt",[], "varfreq",[]);
        end
        
        FilterSummary(end-cycleBackCounter+1).LAfiltered = LAfilt;
        
        clear countArrs altAlleleChar altFrequency filterIndels filterTotalCounts ...
            filterCoverage filterCycleOcc filterAll snp;
    end
    clear c;

    %% Fill SNPSummary
    
    fprintf("\tWrite to summary...\n");

    for i = 1 : length(FilterSummary)
        if ~isempty([FilterSummary(i).filtered.pos])
            SNPSummary(allCount).IndvMutList = vertcat([FilterSummary(i).filtered.pos],string({FilterSummary(i).filtered.ref}),string({FilterSummary(i).filtered.alt}))';
        else 
            SNPSummary(allCount).IndvMutList = [];
        end
        if strcmp(sampleType, "standard")
            SNPSummary(allCount).Sample = sampleCurr+num2str(cycles(i));
        elseif strcmp(sampleType, "customized")
            SNPSummary(allCount).Sample = sampleCurr;
        end
        
        SNPSummary(allCount).LAfiltered = FilterSummary(i).LAfiltered;
        SNPSummary(allCount).FilterSummary = FilterSummary(i);
        allCount = allCount + 1;
    end

end

%% Save data


if saveSNPSummary
    fprintf("Save SNPSummary to '%s' ...\n", savepath + datestr(now, 'yyyymmdd') + '_' + evolExp + "_SNPSummary.mat");
    
    if strcmp(sampleType, "standard")
        save(savepath + datestr(now, 'yyyymmdd') + '_' + evolExp + "_SNPSummary.mat", 'SNPSummary')
    elseif strcmp(sampleType, "customized")
       
            save(savepath + datestr(now, 'yyyymmdd') + "_25comp_SNPSummary.mat", 'SNPSummary')
     
    end
end

disp("Done.");