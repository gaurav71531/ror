function processRiverData(varargin)

% % % declare
global rivPt;
global demPt;
global delXInPdf;
global epsilon;
global lambda;
global flowPowerFactor;

% % % initialize
rivPt = [];
demPt = [];
delXInPdf = 0.05;
epsilon = 1e-9;
if nargin > 0
    lambda = varargin{1};
else
    lambda = 0.01;
end
% lambda = 0;
% flowPowerFactor = efficiency * height drop * gravitational const * (ft^3
% to m^3) / 1000 (for kW)
flowPowerFactor = 0.8 * 8 * 9.8 * 0.3038^3 / 1000;

% % % % Time period for river Data %%%%%%%%%%%
mon = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
dStart = [6,0,0];
dEnd = [6,0,0];
% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % utility functions

% arg: 0 -> extract and save
% arg: 1 -> process riv and dem points
% arg: 2 -> run ror placement algorithm
arg = 2;

switch arg
    case 0
        extractAndSave(dStart, dEnd);
        
    case 1
        processRiverPt(mon{dStart(1)});
        processDemandPt(mon{dStart(1)});
        
    case 2
        % currentData is original without std
        % currentDataNew is with some modifications in pdf mean and std
        % currentDataNew1.mat is with Nebraska
        % currentDataNewTest.mat is with fft values

%         matObj = matfile('currentDataNew1.mat');
%         matObj = matfile('currentDataMar.mat');
        matObj = matfile('currentDataJun.mat');
%         fStr = sprintf('currentData%sWithFFT1.mat', mon{dStart(1)});
%         matObj = matfile(fStr);
        rivPt = matObj.rivPt;
        demPt = matObj.demPt;
        weakRivPtInd = matObj.weakRivPtInd;
        weakRivPtIndNE = matObj.weakRivPtIndNE;
%         neInd = matObj.neInd;

%         weakRivPtInd = [3,6,7,11,14, 16:18, 21, 29:34, 40,45:46, 55, 62:65, 69:83, 86:91, 94, 99, 103, 105, 107, 108, 115, 126, 129, 130, 135, 138, 141:145, 148, 152, 160:168];
%         weakRivPtInd = [155,174:182];
%         weakRivPtIndNE = [140,141];
        % neInd = [140:153];%nebraska 

        rivPt([weakRivPtInd,weakRivPtIndNE]) = [];

%         demPt = demPt(1:2:end);
        progress();
        partitionSpaceForRor(demPt, rivPt, mon{dStart(1)});
        
            
end
% plotRiverCoord();


function partitionSpaceForRor(demPt, rivPt, mon)

tic;
timeBegin = toc;
timeInit = timeBegin;
% global rivPt;
% global demPt;
global lambda;

numDemPt = length(demPt);
% numDemPt = 10;
numRiverPt = length(rivPt);
rorAvailableSet = 1:length(rivPt);
rorOccupiedSet = [];
rorAllocated = 0;
K = 10;
D = cell(K,1);
distCluster = zeros(K,1);
clusterTotDemPdf = struct('x', {}, 'y', {}, 'nzStartInd', {}, 'nzEndInd', {}, 'mean', {}, 'std', {});
% clusterTotDemPdf = struct('x', {}, 'y', {}, 'nzStartInd', {}, 'nzEndInd', {}, 'mean', {}, 'std', {}, 'fft', {});
% clusterTotDemPdf(K) = struct('pdf', []); % allocate empty struct
% distTotal = sum(distCluster);
demPdfTemp = [demPt.pdf];
rivPdfTemp = [rivPt.pdf];
demMeanPlusStdVec = [demPdfTemp.mean] + [demPdfTemp.std];
rivMeanPlusStdVec = [rivPdfTemp.mean] + [rivPdfTemp.std];
% % % LOG
fid = fopen('runLog_full.txt', 'w');

for dInd = 1:numDemPt
%     first allocation
%     if dInd == 107
%         disp('gg');
%     end
    distOfOccupied = zeros(rorAllocated,1);
    demPtAssociationFlag = 0;
    clusterTotDemPdfTemp = struct('x', {}, 'y', {}, 'nzStartInd', {}, 'nzEndInd', {}, 'mean', {}, 'std', {});
%     clusterTotDemPdfTemp = struct('x', {}, 'y', {}, 'nzStartInd', {}, 'nzEndInd', {}, 'mean', {}, 'std', {}, 'fft', {});
    for occRorInd = 1:rorAllocated
%         demClusterInd = [];
%         for dInClusterInd = 1:length(D{occRorInd})
%             demClusterInd = [demClusterInd, D{occRorInd}(dInClusterInd)];
%         end
%         demClusterInd = [demClusterInd,dInd];
        demClusterInd = [D{occRorInd},dInd];
%         if checkValidAssociation(demClusterInd, rivPt(rorOccupiedSet(occRorInd)))
        if checkValidAssociation(demMeanPlusStdVec(demClusterInd), rivMeanPlusStdVec(rorOccupiedSet(occRorInd)))
            [distOfOccupied(occRorInd),clusterTotDemPdfTemp(occRorInd)] = getTotalDistance(demPt(demClusterInd), rivPt(rorOccupiedSet(occRorInd)), 'ind', []);
            demPtAssociationFlag = 1;
        else
            distOfOccupied(occRorInd) = Inf;
        end
    end
    distOfAvailable = zeros(numRiverPt - rorAllocated,1);
    for rorAvailInd = 1:length(rorAvailableSet)
        if checkValidAssociation(demMeanPlusStdVec(dInd), rivMeanPlusStdVec(rorAvailableSet(rorAvailInd)))
%         if checkValidAssociation(dInd, rivPt(rorAvailableSet(rorAvailInd)))
            distOfAvailable(rorAvailInd) = getTotalDistance(demPt(dInd), rivPt(rorAvailableSet(rorAvailInd)), 'ind', []);
            demPtAssociationFlag = 1;
        else
            distOfAvailable(rorAvailInd) = Inf;
        end
    end
    
%     first assignment
    if ~demPtAssociationFlag
        fprintf('demPt-%d cannot be assigned to any ROR, the program will terminate\n', dInd);
        return;
    end
    [dMinAvail, indMinAvail] = min(distOfAvailable);
    [~, indMinOcc] = min(distOfOccupied - distCluster(1:rorAllocated));
    distClusterWAvail = distCluster;
    distClusterWOcc = distCluster;
    
    if rorAllocated < K
        rorAllocatedTemp = rorAllocated + 1;
        distClusterWAvail(rorAllocatedTemp) = dMinAvail;
    else
        distClusterWAvail(end) = Inf; % top stop new ror pickup
    end
    
    if isempty(distOfAvailable), distClusterWAvail(1) = Inf;end
    if isempty(distOfOccupied)
        distClusterWOcc(1) = Inf;
    else
        distClusterWOcc(indMinOcc) = distOfOccupied(indMinOcc);
    end
    
    if ~(sum(distClusterWAvail) == Inf && sum(distClusterWOcc) == Inf)
        
        if sum(distClusterWAvail) < sum(distClusterWOcc)
    %         assign to new ror site
            rorAllocated = rorAllocatedTemp;
            D{rorAllocated} = [D{rorAllocated}, dInd];
            rorOccupiedSet = [rorOccupiedSet, rorAvailableSet(indMinAvail)];
    % % %         LOG
            fprintf(fid, 'demPt-%d assigned to new ROR-%d\n', dInd, rorAvailableSet(indMinAvail));
            rorAvailableSet(indMinAvail) = [];        
            distCluster = distClusterWAvail;
            clusterTotDemPdf(rorAllocated) = demPt(dInd).pdf;
        else
    %         assign to existing ones
            D{indMinOcc} = [D{indMinOcc}, dInd];
            distCluster = distClusterWOcc;
            clusterTotDemPdf(indMinOcc) = clusterTotDemPdfTemp(indMinOcc);
            % % % LOG
            fprintf(fid, 'demPt-%d assigned to existing ROR-%d\n', dInd, rorOccupiedSet(indMinOcc));
        end

        distTotal = sum(distCluster);
        %   intra- group shuffle
        
        for rorInd = 1:rorAllocated
            if length(D{rorInd}) == 1, continue;end
            distRorAvailable = zeros(numRiverPt - rorAllocated, 1);
            for rorAvailableInd = 1:length(rorAvailableSet)
                if checkValidAssociation(demMeanPlusStdVec(D{rorInd}), rivMeanPlusStdVec(rorAvailableSet(rorAvailableInd)))
%                 if checkValidAssociation(D{rorInd}, rivPt(rorAvailableSet(rorAvailableInd)))
%                     if length(clusterTotDemPdf(rorInd).nzStartInd) >1
%                     end
                    [distRorAvailable(rorAvailableInd), ~] = getTotalDistance(demPt(D{rorInd}), rivPt(rorAvailableSet(rorAvailableInd)), 'pdf', clusterTotDemPdf(rorInd));
                else
                    distRorAvailable(rorAvailableInd) = Inf;
                end     
            end
    %         distRorOccupied = zeros(rorAllocated, 1);
    %         for rorOccupiedInd = 1:rorAllocated
    %             if rorOccupiedSet(rorInd) == rorOccupiedSet(rorOccupiedInd),continue;end
    %             distRorOccupied(rorOccupiedInd) = getTotalDistance([D{rorInd},D{rorOccupiedInd}], rivPt(rorAvailableSet(rorAvailableInd)), 'ind');
    %         end
            [dMinAvail, indMinAvail] = min(distRorAvailable);
    %         [dMinOcc, indMinOcc] = min(distRorOccupied);
            distClusterOcc = ones(length(distCluster),1) * Inf;
            distClusterAvail = distCluster;
    %         distClusterOcc(rorInd) = 0;
    %         distClusterOcc(indMinAvail) = dMinOcc;

            distClusterAvail(rorInd) = dMinAvail;
            sumDistClusterAvail = sum(distClusterAvail);
            sumDistClusterOcc = sum(distClusterOcc);
            dMin = min(sumDistClusterAvail, sumDistClusterOcc);
            if dMin < distTotal
                if sumDistClusterAvail < sumDistClusterOcc
                    indReturn = rorOccupiedSet(rorInd);
                    fprintf(fid, 'group with ror-%d moved to site-%d\n', indReturn, rorAvailableSet(indMinAvail));
                    rorOccupiedSet(rorInd) = rorAvailableSet(indMinAvail);
                    rorAvailableSet(indMinAvail) = indReturn;
                    distCluster(rorInd) = distRorAvailable(indMinAvail);
                else
    %                 not used right now
                    rorAvailableSet = [rorAvailableSet, rorOccupiedSet(rorInd)];
                    rorOccupiedSet(rorInd) = [];
                    D{indMinOcc} = [D{indMinOcc}, D{rorInd}];
                    D{rorInd} = [];
                end
            end
        end
    else
%         no new ROR can be used && existing ROR cannot afford this point
%         This would require entire group to look for new ROR-site
        reshuffleFlag = 0;
        distClusterReshuffleCurrent = Inf;
        clusterTotDemPdfTemp = clusterTotDemPdf;
%         if dInd==76
%         end
        for rorInd = 1:rorAllocated
            DTemp = D;
            DTemp{rorInd} = [DTemp{rorInd}, dInd];
            
            distRorAvailable = zeros(numRiverPt - rorAllocated, 1);
            for rorAvailableInd = 1:length(rorAvailableSet)
                if checkValidAssociation(demMeanPlusStdVec(DTemp{rorInd}), rivMeanPlusStdVec(rorAvailableSet(rorAvailableInd)))
%                 if checkValidAssociation(DTemp{rorInd}, rivPt(rorAvailableSet(rorAvailableInd)))
                    [distRorAvailable(rorAvailableInd), clusterTotDemPdfTemp(rorInd)] = getTotalDistance(demPt(DTemp{rorInd}), rivPt(rorAvailableSet(rorAvailableInd)), 'ind', clusterTotDemPdf(rorInd));
                else
                    distRorAvailable(rorAvailableInd) = Inf;
                end     
            end
            [dMinAvail, indMinAvail] = min(distRorAvailable);
            distClusterAvail = distCluster;
            
            distClusterAvail(rorInd) = dMinAvail;
            sumDistClusterAvail = sum(distClusterAvail);
            
            if sumDistClusterAvail < distClusterReshuffleCurrent
                distClusterMin = distClusterAvail;
                distClusterReshuffleCurrent = sumDistClusterAvail;
                rorAvailableSetMin = rorAvailableSet;
                rorOccupiedSetMin = rorOccupiedSet;
                indReturn = rorOccupiedSet(rorInd);
                rorOccupiedSetMin(rorInd) = rorAvailableSet(indMinAvail);
                rorAvailableSetMin(indMinAvail) = indReturn;
                DMin = DTemp;

                indReshuffleMin = indReturn;
                reshuffleFlag = 1;
                newSiteNum = rorAvailableSet(indMinAvail);
                
                clusterTotDemPdfTempMin = clusterTotDemPdfTemp;
            end
        end
        if reshuffleFlag
            distCluster = distClusterMin;
            rorAvailableSet = rorAvailableSetMin;
            rorOccupiedSet = rorOccupiedSetMin;
            D = DMin;
            clusterTotDemPdf = clusterTotDemPdfTempMin;
            fprintf(fid, 'group with ROR-%d moved to site-%d due to insufficiency\n', indReshuffleMin, newSiteNum);
        else
            fprintf('New demPt-%d cannot be associated anywhere and anyhow by the program, will termonate\n', dInd);
        end        
    end
    
    
%     coordinate descent (reshuffling)

    distTotal = sum(distCluster);

%   inter-group shuffle
    for dInOccRorInd = 1:dInd
%         try
        occRorInd = locateDemPtInRorCluster(dInOccRorInd, D);
%         catch
%         end
        
        DTempRef = D;
        distClusterRef = distCluster;
        clusterTotDemPdfRef = clusterTotDemPdf;
        DTempRef{occRorInd}(DTempRef{occRorInd} == dInOccRorInd) = [];
%         demClusterInd = [];
%         for i = 1:length(DTempRef{occRorInd})
%             demClusterInd = [demClusterInd, DTempRef{occRorInd}(i)];
%         end
        demClusterInd = DTempRef{occRorInd};
        if checkValidAssociation(demMeanPlusStdVec(demClusterInd), rivMeanPlusStdVec(rorOccupiedSet(occRorInd)))
%         if checkValidAssociation(demClusterInd, rivPt(rorOccupiedSet(occRorInd)))
            [distClusterRef(occRorInd), clusterTotDemPdfRef(occRorInd)] = getTotalDistance(demPt(demClusterInd), rivPt(rorOccupiedSet(occRorInd)), 'ind', []);
        else
            disp('some error here');
            distClusterRef(occRorInd) = Inf;
        end
        distTotalCurrent = distTotal;
        reshuffleFlag = 0;
        for occRorIndNext = 1:rorAllocated
            if (occRorIndNext == occRorInd), continue;end
            DTemp = DTempRef;
            DTemp{occRorIndNext} = [DTemp{occRorIndNext}, dInOccRorInd];
%             demClusterInd = [];
%             for i = 1:length(DTemp{occRorIndNext})
%                 demClusterInd = [demClusterInd, DTemp{occRorIndNext}(i)];
%             end
            demClusterInd = DTemp{occRorIndNext};
            distClusterTemp = distClusterRef;
            clusterTotDemPdfTemp = clusterTotDemPdfRef;
            if checkValidAssociation(demMeanPlusStdVec(demClusterInd), rivMeanPlusStdVec(rorOccupiedSet(occRorIndNext)))
%             if checkValidAssociation(demClusterInd, rivPt(rorOccupiedSet(occRorIndNext)))
                [distClusterTemp(occRorIndNext), clusterTotDemPdfTemp(occRorIndNext)] = getTotalDistance(demPt(demClusterInd), rivPt(rorOccupiedSet(occRorIndNext)), 'ind', []);
            else
                distClusterTemp(occRorIndNext) = Inf;
            end
            if sum(distClusterTemp) == 0
                fprintf('zero sum for distance vector\n');
            end
            if sum(distClusterTemp) < distTotalCurrent
                distTotalCurrent = sum(distClusterTemp);
                reshuffleNextInd = occRorIndNext;
                reshuffleFlag = 1;
                distClusterMin = distClusterTemp;
                clusterTotDemPdfMin = clusterTotDemPdfTemp;
            end
        end
        if reshuffleFlag
            distCluster = distClusterMin;
            clusterTotDemPdf = clusterTotDemPdfMin;
            distTotal = sum(distCluster);
            DTemp = D;
            DTemp{occRorInd}(DTemp{occRorInd} == dInOccRorInd) = [];
            DTemp{reshuffleNextInd} = [DTemp{reshuffleNextInd}, dInOccRorInd];
            freeInd = zeros(1,rorAllocated);
            % % % LOG
            fprintf(fid, 'demPt-%d reshuffled to ROR-%d\n', dInOccRorInd, rorOccupiedSet(reshuffleNextInd));
            indUse = 1;
            for i = 1:rorAllocated
                if isempty(DTemp{i})
                    freeInd(indUse) = i;
                    indUse = indUse + 1;
                    % % % LOG
                    fprintf(fid, 'ROR-%d is left empty\n', rorOccupiedSet(i));
                end
            end
            freeInd(indUse:end) = [];
            if ~isempty(freeInd)
                DTemp(freeInd) = [];
                DTempNext = cell(K,1);
                DTempNext(1:length(DTemp)) = DTemp;
                DTemp = DTempNext;
                rivPtFree = rorOccupiedSet(freeInd);
                rorOccupiedSet(freeInd) = [];
                rorAvailableSet = [rorAvailableSet, rivPtFree];
                rorAllocated = rorAllocated - length(freeInd);

                distCluster(freeInd) = [];
                distClusterNext = zeros(K,1);
                distClusterNext(1:length(distCluster)) = distCluster;
                distCluster = distClusterNext;
                distTotal = sum(distCluster);
                
                clusterTotDemPdf(freeInd) = [];
            end
            D = DTemp;
        end
    end 
    progress(dInd / numDemPt);
    timeFin = toc;
    if timeFin - timeInit > 300
        fprintf('simulation completed = %d %%, time elapsed = %f\n', floor(dInd * 100 / numDemPt), timeFin - timeBegin);
        timeInit = timeFin;
    end
end
% % LOG
fclose(fid);
fprintf('simulated completed = 100%%, time elepased = %f\n', toc-timeBegin);

assignin('base', 'D', D);
assignin('base', 'rorOccupiedSet', rorOccupiedSet);
assignin('base', 'distCluster', distCluster);
assignin('base', 'rivPt', rivPt);
assignin('base', 'demPt', demPt);
assignin('base', 'lambda', lambda);
fName = sprintf('/OutputData/output_lambda_%4.3f_numDemPt_%d_%sTemp.mat', lambda, length(demPt), mon);
save([pwd fName], 'D', 'distCluster', 'rorOccupiedSet', 'rivPt', 'demPt', 'lambda');
% plotOutputData(D, rorOccupiedSet, demPt, lambda);


% function [rorAvailableSet,rorOccupiedSet,distCluster,newSiteNum] = performGroupReshuffle(rorAllocated, D, numRiverPt, rorAvailableSet, rorOccupiedSet, clusterTotDemPdf, distTotal, distCluster, fid, type)
% 
% % global demPt;
% global rivPt;
% 
% newSiteNum = 999999; % dummy large value
% 
% for rorInd = 1:rorAllocated
%     if length(D{rorInd}) == 1, continue;end
%     distRorAvailable = zeros(numRiverPt - rorAllocated, 1);
%     for rorAvailableInd = 1:length(rorAvailableSet)
%         if checkValidAssociation(D{rorInd}, rivPt(rorAvailableSet(rorAvailableInd)))
%             [distRorAvailable(rorAvailableInd), ~] = getTotalDistance(D{rorInd}, rivPt(rorAvailableSet(rorAvailableInd)), 'pdf', clusterTotDemPdf(rorInd));
%         else
%             distRorAvailable(rorAvailableInd) = Inf;
%         end     
%     end
% %         distRorOccupied = zeros(rorAllocated, 1);
% %         for rorOccupiedInd = 1:rorAllocated
% %             if rorOccupiedSet(rorInd) == rorOccupiedSet(rorOccupiedInd),continue;end
% %             distRorOccupied(rorOccupiedInd) = getTotalDistance([D{rorInd},D{rorOccupiedInd}], rivPt(rorAvailableSet(rorAvailableInd)), 'ind');
% %         end
%     [dMinAvail, indMinAvail] = min(distRorAvailable);
% %         [dMinOcc, indMinOcc] = min(distRorOccupied);
%     distClusterOcc = ones(length(distCluster),1) * Inf;
%     distClusterAvail = distCluster;
% %         distClusterOcc(rorInd) = 0;
% %         distClusterOcc(indMinAvail) = dMinOcc;
% 
%     distClusterAvail(rorInd) = dMinAvail;
%     sumDistClusterAvail = sum(distClusterAvail);
%     sumDistClusterOcc = sum(distClusterOcc);
%     dMin = min(sumDistClusterAvail, sumDistClusterOcc);
%     if dMin < distTotal
%         if sumDistClusterAvail < sumDistClusterOcc
%             indReturn = rorOccupiedSet(rorInd);
%             if strcmp(type, 'run')
%                 fprintf(fid, 'group with ror-%d moved to site-%d\n', indReturn, rorAvailableSet(indMinAvail));
%             end
%             newSiteNum = rorAvailableSet(indMinAvail);
%             rorOccupiedSet(rorInd) = rorAvailableSet(indMinAvail);
%             rorAvailableSet(indMinAvail) = indReturn;
%             distCluster(rorInd) = distRorAvailable(indMinAvail);
%         else
% %                 not used right now
%             fprintf('should not reach here\n');
% %             rorAvailableSet = [rorAvailableSet, rorOccupiedSet(rorInd)];
% %             rorOccupiedSet(rorInd) = [];
% %             D{indMinOcc} = [D{indMinOcc}, D{rorInd}];
% %             D{rorInd} = [];
%         end
%     end
% end


function [val, demPdf] = getTotalDistance(demClusterPt, rivPt, type, pdfUse)

% global demPt;
global lambda;

demPdf = struct('x', [], 'y', [], 'nzStartInd', [], 'nzEndInd', [], 'mean', [], 'std', []);
% demPdf = struct('x', [], 'y', [], 'nzStartInd', [], 'nzEndInd', [], 'mean', [], 'std', [], 'fft', []);

if isempty(demClusterPt)
    val = 0;
    return;
end
euclideanDist = 0;
for i = 1:length(demClusterPt)
    euclideanDist = euclideanDist + (demClusterPt(i).xPos - rivPt.xPos)^2 ...
        + (demClusterPt(i).yPos - rivPt.yPos)^2;
end

if strcmp(type, 'ind') 
    if length(demClusterPt) == 1
        demPdf = demClusterPt(1).pdf;
    else
        pdf1 = demClusterPt(1).pdf;
        for i = 2:length(demClusterPt)
            pdf2 = demClusterPt(i).pdf;
            demPdf = convolve(pdf1, pdf2);
            pdf1 = demPdf;
        end
%         demPdf = computerDemPdfParallel(demClusterPt);
    end
elseif strcmp(type, 'pdf')
    demPdf = pdfUse;
end
% demPdf.fft = [];
klDist = getKLDist(demPdf, rivPt.pdf);

val = klDist + lambda* euclideanDist;


function demPdf = computerDemPdfParallel(demPt)

if length(demPt) == 1
    demPdf = demPt(1).pdf;
    return
end

% demPtNew = struct('xPos', {}, 'yPos', {}, 'pdf', {});
% demPtNew(length(demPt)) = struct('xPos', [], 'yPos', [], 'pdf', []);
if rem(length(demPt),2)
    lenBy2 = (length(demPt)-1)/2;
    demPt1 = demPt(1:lenBy2);
    demPt2 = demPt(lenBy2+1:lenBy2*2);
    demPtNew = struct('xPos', cell(1,length(lenBy2+1)), 'yPos', cell(1,length(lenBy2+1)), 'pdf', cell(1,length(lenBy2+1)));
    for i = 1:lenBy2
        demPtTemp = struct('xPos', [], 'yPos', [], 'pdf', convolve(demPt1(i).pdf, demPt2(i).pdf));
        demPtNew(i) = demPtTemp;
%         demPtNew(i).pdf = convolve(demPt1(i).pdf, demPt2(i).pdf);
    end
    demPtNew(lenBy2+1).pdf = demPt(end).pdf;
else
    lenBy2 = length(demPt)/2;
    demPt1 = demPt(1:lenBy2);
    demPt2 = demPt(lenBy2+1:end);
    demPtNew = struct('xPos', cell(1,length(lenBy2)), 'yPos', cell(1,length(lenBy2)), 'pdf', cell(1,length(lenBy2)));
    for i = 1:lenBy2
        demPtTemp = struct('xPos', [], 'yPos', [], 'pdf', convolve(demPt1(i).pdf, demPt2(i).pdf));
        demPtNew(i) = demPtTemp;
%         demPtNew(i).pdf = convolve(demPt1(i).pdf, demPt2(i).pdf);
    end
end

demPdf = computerDemPdfParallel(demPtNew);


function val = getKLDist(p, q)


global delXInPdf;
global epsilon;

% val = 0;
% 
% for pLocalInd = 1:length(p.y)
%     actualInd = p.nzStartInd + pLocalInd - 1;
% %     try
%     if actualInd >= q.nzStartInd && actualInd <= q.nzEndInd
%         qLocalInd = actualInd - q.nzStartInd + 1;
%     end
% %     catch
% %     end
%     if actualInd < q.nzStartInd || actualInd > q.nzEndInd
%         qValUse = epsilon;
%     else % in between closed interval 
%         if q.y(qLocalInd) < epsilon
%             qValUse = epsilon;
%         else
%             qValUse = q.y(qLocalInd);
%         end
%     end
%     val = val + p.y(pLocalInd) * log2(p.y(pLocalInd) / qValUse) * delXInPdf;
% end
if p.nzStartInd- q.nzStartInd + 1 <= 0
    qLocalStartInd = 1;
    pLocalOffsetInd = -p.nzStartInd + q.nzStartInd + 1;
else
    qLocalStartInd = p.nzStartInd - q.nzStartInd + 1;
    pLocalOffsetInd = 1;
end

qLocalStopInd = min(p.nzStartInd + length(p.y) - 1, q.nzEndInd) - q.nzStartInd + 1;

qUse = q.y(qLocalStartInd:qLocalStopInd);
qUse(qUse < epsilon) = epsilon;
qTemp = ones(1,length(p.y))*epsilon;
qTemp(pLocalOffsetInd:pLocalOffsetInd+length(qUse)-1) = qUse;

val = sum(p.y .* log2(p.y ./ qTemp)) * delXInPdf;
% 
if isnan(val), val = Inf;end


function out = convolve(fcn1, fcn2)

global delXInPdf;
% needs struct with x, y and non-zero interval [nzStartInd, nzEndInd]
% out = struct('x', 0, 'y', 0, 'nzStartInd', 0, 'nzEndInd', 0, 'mean', 0, 'std', 0);
out.x = 0;
% out.y = 0;
% x1 and x2 must have same starting point and same step-size
% if (fcn1.x(1) ~= fcn2.x(1)), 
%     fprintf('error in input x-axis format\n');
%     return;
% end
% l1 = length(fcn1.y);
% l2 = length(fcn2.y);
% yOld = conv(fcn1.y, fcn2.y) * delX;

% out.x = [0:delX:delX*((l1+l2-1)-1)] + fcn1.x(1);

% y = conv(fcn1.y(fcn1.nzStartInd:fcn1.nzEndInd), fcn2.y(fcn2.nzStartInd:fcn2.nzEndInd));
numZerosBegin = fcn1.nzStartInd + fcn2.nzStartInd - 2;
% numZerosEnd = l1 - fcn1.nzEndInd + l2 - fcn2.nzEndInd; 
% out.y = [zeros(1,numZerosBegin), y, zeros(1, numZerosEnd)] * delX; % zero padding on beginning and end

% [sTemp, eTemp] = getNZInterval(y);
% out.nzStartInd = numZerosBegin + sTemp;
% out.nzEndInd = numZerosBegin + eTemp;
y = conv(fcn1.y, fcn2.y);
% y = y * delXInPdf;
y = y * 0.05;
[s, e] = getNZInterval(y);
out.nzStartInd = numZerosBegin + s; 
% length(out.nzStartInd)
% if length(out.nzStartInd) > 1
% end
out.y = y(s:e);
out.nzEndInd = out.nzStartInd + length(out.y) - 1;
out.mean = fcn1.mean + fcn2.mean;
out.std = sqrt(fcn1.std^2 + fcn2.std^2);

% if (length(out.y) ~= l1 + l2 -1), 
%     error('error in conv length');
% end

function occRorInd = locateDemPtInRorCluster(dInd, D)


for i = 1:length(D)
    for j = 1:length(D{i})
        if dInd == D{i}(j)
            occRorInd = i;
            return;
        end
    end
end


% function out = checkValidAssociation(demClusterInd, rivPt)
function out = checkValidAssociation(dem, gen)

% global demPt;
out = 0;

% sumMeanStd = 0;
% for i = 1:length(demClusterInd)
%     sumMeanStd = sumMeanStd + demPt(demClusterInd(i)).pdf.mean + demPt(demClusterInd(i)).pdf.std;
% end
if sum(dem) < gen
    out = 1;
end

% if sumMeanStd < rivPt.pdf.mean + rivPt.pdf.std
%     out = 1;
% end


function processRiverPt(month)


global rivPt;
global delXInPdf;

folderStr = sprintf('/Data/sitePowerPdf%s/', month);
% dataFile = {'latLon_Minnesota.mat', 'latLon_Iowa.mat', 'latLon_Illinois.mat', 'latLon_Missouri.mat', 'latLon_Wisconsin.mat'};
% dataFile = {'latLon_Minnesota.mat'};
% colorVec = {'b', 'r', ', 'm', 'c', 'k'};
fileName = dir([pwd folderStr]);

% legendStr = cell(1,length(dataFile));
X = zeros(length(fileName), 1);
Y = zeros(length(fileName), 1);
meanVec = zeros(length(fileName), 1);
vecInd = 1;
siteIDVec = cell(length(fileName),1);
for fileInd = 1:length(fileName)
    if length(fileName(fileInd).name) < 5, continue;end
    if ~strcmp(fileName(fileInd).name(end-2:end), 'mat'), continue;end
    
    matObj = matfile([pwd folderStr, fileName(fileInd).name]);
%     pdfVal = matObj.f;
%     flowVal = matobj.xi;
    latVal = matObj.latVal;
    lonVal = matObj.lonVal;
    meanPowGen = matObj.meanPowGen;
    varPowGen = matObj.varPowGen;
    
    yUse = matObj.f;
    xUse = matObj.xi;
    
    nzOffset = floor(xUse(1)/delXInPdf);
    [~, e] = getNZInterval(yUse);
    if isempty(e), e = 0;end
    pdf(vecInd).nzStartInd = nzOffset + 1;
    pdf(vecInd).nzEndInd = nzOffset + e;
    pdf(vecInd).x = xUse;
    pdf(vecInd).y = yUse(1:e);
    pdf(vecInd).mean = meanPowGen;
    pdf(vecInd).std = sqrt(varPowGen);
    
%     figure;
%     plot(pdf(vecInd).x, pdf(vecInd).y);grid;
%     if meanFlow > 1e4,continue;end
%     varFlow = matObj.varFlow;
%     dczone = utmzone(latVal, lonVal);
    utmstruct = defaultm('utm');
%     utmstruct.zone = dczone;
    utmstruct.zone = '15T'; % plot multi-zone data in single zone for map consistency
    utmstruct.geoid = wgs84Ellipsoid('meters');
    
    utmstruct = defaultm(utmstruct);
    [X(vecInd),Y(vecInd)] = mfwdtran(utmstruct,latVal,lonVal);
    meanVec(vecInd) = meanPowGen;
    siteIDVec{vecInd} = fileName(fileInd).name(1:end-4);
    vecInd = vecInd + 1;
%     hold on;
%     scatter(x, y, colorVec{1}, 'filled');
%      scatter bug 1283854  %
%     scattStr = sprintf('s%d.MarkerEdgeColor = s%d.CData;', i,i);
%     eval(scattStr);
    
%     c = cellstr(num2str(fileInd));

%     text(x, y, c,'color', colorVec{2},'fontsize', 15);
end
X(vecInd:end) = [];
Y(vecInd:end) = [];
meanVec(vecInd:end) = [];
siteIDVec(vecInd:end) = [];

X = X / 1e5;  % normalizing by 100 km
Y = Y / 1e5;

rivPt = struct('xPos', num2cell(X), 'yPos', num2cell(Y), ...
    'meanPowGen', num2cell(meanVec), 'siteId', siteIDVec);
for i = 1:vecInd-1
    rivPt(i).pdf = pdf(i);
end
% rivPt.xPos = X;
% rivPt.yPos = Y;
% rivPt.meanFlw = meanVec;
% rivPt.siteId = siteIDVec;
% assignin('base', 'siteIDVec', siteIDVec);
% assignin('base', 'meanVec', meanVec);


function processDemandPt(month)


global demPt;
global delXInPdf;

fStr = sprintf('/Data/demand/demandProfile%s.mat', month);
matObj = matfile([pwd fStr]);
demPtTemp = matObj.demPt;

x = 0:delXInPdf:3;
for i = 1:length(demPtTemp.x)
    demPt(i).xPos = demPtTemp.x(i);
    demPt(i).yPos = demPtTemp.y(i);
    
    pdfObj.x = x;
    yFcn = getGaussMixModFn([demPtTemp.mu1(i), demPtTemp.mu2(i)], [demPtTemp.sigSq1(i), demPtTemp.sigSq2(i)], demPtTemp.alpha(i));
    yUse = yFcn(pdfObj.x);
    [pdfObj.nzStartInd, pdfObj.nzEndInd] = getNZInterval(yUse);
    if length(pdfObj.nzStartInd) > 1
    end
    pdfObj.y = yUse(pdfObj.nzStartInd:pdfObj.nzEndInd);
    demPt(i).pdf = pdfObj;
    mean = sum(x .* yUse) * delXInPdf;
    var = sum((x - mean).^2 .* yUse) * delXInPdf;
    demPt(i).pdf.mean = mean;
    demPt(i).pdf.std = sqrt(var);
end

% dCoord.yPos = demPt.y;

% dCoord.mu1 = demPt.mu1;
% dCoord.mu2 = demPt.mu2;

% dCoord.sigSq1 = demPt.sigSq1;
% dCoord.sigSq2 = demPt.sigSq2;

% dCoord.alpha = demPt.alpha;
% N = 10;
% for i = 1:N
%     demPt(i).x = x;
%     yFcn = getGaussMixModFn([dCoord.mu1(i), dCoord.mu2(i)], [dCoord.sigSq1(i), dCoord.sigSq2(i)], dCoord.alpha(i));
%     demPt(i).y = yFcn(x);
%     [demPt(i).nzStartInd, demPt(i).nzEndInd] = getNZInterval(demPt(i).y);
% end


function y = getGaussMixModFn(mu, sigSq, alpha)


y = @(x) alpha * 1/sqrt(2*pi*sigSq(1))*exp(-(x-mu(1)).^2/2/sigSq(1)) ...
    + (1-alpha) * 1/sqrt(2*pi*sigSq(2))*exp(-(x-mu(2)).^2/2/sigSq(2));


function [s, e] = getNZInterval(y)

% global epsilon;
epsilon = 1e-9;

nz = y >= epsilon;
s = find(nz, 1, 'first');
e = find(nz, 1, 'last');


function plotRiverCoord()


global rivPt;
global demPt;
% dataFile = {'latLon_Minnesota.mat', 'latLon_Iowa.mat', 'latLon_Illinois.mat', 'latLon_Missouri.mat', 'latLon_Wisconsin.mat'};
% dataFile = {'latLon_Minnesota.mat'};
colorVec = {'b', 'r', 'g', 'm', 'c', 'k'};

figure;
colormap jet;
scatter([rivPt.xPos], [rivPt.yPos], 35, [rivPt.meanPowGen], 'filled');
c = cellstr(num2str([1:length(rivPt)]'));
text([rivPt.xPos], [rivPt.yPos], c, 'color', 'k','fontsize', 12);
colorbar;
xlabel('x (\times 100 km)');
ylabel('y (\times 100 km)');
grid;

% adding demand pts
hold on;

scatter([demPt.xPos], [demPt.yPos], 10, 'r');
c = cellstr(num2str([1:length(demPt)]'));
hold on;
text([demPt.xPos], [demPt.yPos], c, 'color', 'r','fontsize', 10);


matObj = matfile([pwd '/Data/demand/demandVars.mat']);
xAxesVal = matObj.xAxesVal;
yAxesVal = matObj.yAxesVal;
axis([xAxesVal yAxesVal]);


function plotPdf(demCluster)

global delXInPdf;

if length(demCluster) > 1
    pdf1 = demCluster(1).pdf;
    for i = 2:length(demCluster)
        pdf2 = demCluster(i).pdf;
        demPdf = convolve(pdf1, pdf2);
        pdf1 = demPdf;
    end
    
else
    pdf1 = demCluster.pdf;
end

figure;
x = delXInPdf*(pdf1.nzStartInd-1)+[0:delXInPdf:delXInPdf*(length(pdf1.y)-1)];
plot(x, pdf1.y);grid;


function extractAndSave(dStart, dEnd)


global delXInPdf;
global flowPowerFactor;
% fileName = 'test39.csv';
fileName = dir([pwd '/Data/siteFlow']);
% fileName.name = 'flow_siteID_05287890.csv';
% dStart = [3,0,0];
% dEnd = [3,0,0];
% dStart = [6,0,0];
% dEnd = [6,0,0];
mon = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
% % % % % 
folderStr = sprintf('/Data/sitePowerPdf%s/', mon{dStart(1)});
mkdir([pwd folderStr]);
% % % % % % % % % % % % % % 
% %first week of the year
% dStart = [1,1,0];
% dEnd = [1,7,0];
%month of march

% % % % % % % % % % % % % % 
% siteIDNE = {'06805600','06806500','06807000','06810070','06811500','06813500','06814500',',06815000','06466010','06466400','06466500','06467500','06478522','06478523','06478526','06486000','06600900','06601000','06601200','06609100','06610000','06610705','06610710', '06610720'};

for fileInd = 1:length(fileName)
    fileNameUse = fileName(fileInd).name;
    if length(fileNameUse) < 4, continue;end
    if ~strcmp(fileNameUse(1:4), 'flow'), continue;end
    siteID = fileNameUse(13:end-4);
%     count = 0;
%     for iGG = 1:length(siteIDNE)
%         if strcmp(siteID, siteIDNE{iGG})
%             count = count + 1;
%         end
%     end
%     if count < 1,continue;end
    fprintf('data read begin for siteID=%s ... ', siteID); 
    fid = fopen([pwd, '/Data/siteFlow/', fileName(fileInd).name]);
    dataRead = textscan(fid, '%s %s %*s %s', 'CommentStyle', '#'); % skipping third column, time zone
    fclose(fid);
    
    data.date = dataRead{1};
    data.time = dataRead{2};
    data.flow = dataRead{3};
    
    numData = size(data.date,1);
    
    flowUse = zeros(numData,1);
    flowInd = 1;
    for indData = 1:numData
        if isempty(data.date{indData}) || isempty(data.time{indData}) || isempty(data.flow{indData})
            continue;
        end
        if IsDateInBetween(data.date{indData},dStart,dEnd, 'month')
            valTemp = str2double(data.flow{indData}(1:end-1));
            if isnan(valTemp), continue;end
            flowUse(flowInd) = valTemp;
            flowInd = flowInd + 1;
        end
    end
    flowUse(flowInd:end) = [];
    if isempty(flowUse), continue;end
    powGen = flowPowerFactor * flowUse;
    
    [~, xiTemp] = ksdensity(powGen);
    
    xi = floor(min(xiTemp)/delXInPdf) * delXInPdf : delXInPdf : ceil(max(xiTemp)/delXInPdf) * delXInPdf;
    [f,~] = ksdensity(powGen, xi);    
    
%     xiGG = floor(min(powGen)/delXInPdf) * delXInPdf : delXInPdf : ceil(max(powGen)/delXInPdf) * delXInPdf;
%     [fGG,~] = ksdensity(powGen, xiGG);    
%     figure;
%     plot(xiTemp, fTemp);
%     hold on;
%     plot(xi, f);
%     hold on;
%     plot(xiGG, fGG);
    
    meanPowGen = mean(powGen);
    varPowGen = var(powGen);
    fid = fopen([pwd, '/Data/siteLatLon/siteID_Loc_', siteID, '.txt']);
    dataLL = textscan(fid, '%f %f');
    latVal = dataLL{1};
    lonVal = dataLL{2};
    save([pwd, folderStr, siteID, '.mat'], 'f', 'xi', 'meanPowGen', 'varPowGen', 'latVal', 'lonVal');
    fprintf('data saved for siteID = %s\n', siteID);
end


function out = IsDateInBetween(date, dStart, dEnd, typ)

% check if the day is in between including boundaries
dateUse = textscan(date, '%d', 'delimiter', '/');

out = 0;
if strcmp(typ, 'week')
    if (dateUse{1}(1) >= dStart(1) && dateUse{1}(1) <= dEnd(1)) %month
        if (dateUse{1}(2) >= dStart(2) && dateUse{1}(2) <= dEnd(2)) %day
            out = 1;
        end
    end
elseif strcmp(typ, 'month')
    
    if (dateUse{1}(1) == dStart(1)) %month
        out = 1;
    end
end

% function plotOutputData(varargin)
% 
% global demPt;
% global rivPt;
% global lambda;
% 
% colorVec = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1], [0.7 0.1 0.9], [0 0 0], [1, 0.5, 0.5], [0.5, 0.5, 0.5]};
% % colorVec = {'b', 'r', 'g', 'm', 'c', 'k', 'y', 'w'};
% 
% matObj = matfile('currentData.mat');
% rivPt = matObj.rivPt;
% % demPt = matObj.demPt;
% % demPt = demPt(1:2:end);
% 
% if nargin < 2
%     matObj = matfile([pwd '/OutputData/output_lambda_0.000_numDemPt_576.mat']);
%     D = matObj.D;
%     rorOccupiedSet = matObj.rorOccupiedSet;
%     demPt = matObj.demPt;
% else
%     D = varargin{1};
%     rorOccupiedSet = varargin{2};
%     demPt = varargin{3};
% end
% D = D(1:length(rorOccupiedSet));
%     
% 
% % reorder based on cardinality
% cardSet = cellfun(@(x) length(x), D);
% [~, ind] = sort(cardSet, 'descend');
% D = D(ind);
% rorOccupiedSet = rorOccupiedSet(ind);
% 
% % figure;
% % colormap jet;
% % scatter([rivPt.xPos], [rivPt.yPos], 35, [255,214,211]/255, 'filled');
% % hold on;
% % % scatter([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], 50, 'r', 'filled');
% % c = cellstr(num2str(rorOccupiedSet'));
% % text([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], c, 'color', 'k','fontsize', 12);
% % % colorbar;
% % xlabel('x (\times 100 km)');
% % ylabel('y (\times 100 km)');
% % grid;
% % 
% % % adding demand pts
% % hold on;
% % 
% % for i = 1:length(rorOccupiedSet)
% %     scatter(rivPt(rorOccupiedSet(i)).xPos, rivPt(rorOccupiedSet(i)).yPos, 70, colorVec{i}, 'filled');
% %     hold on;
% %     scatter([demPt(D{i}).xPos], [demPt(D{i}).yPos], 10, colorVec{i});
% %     hold on;
% % end
% % hold off;
% % 
% % matObj = matfile([pwd '/Data/demand/demandVars.mat']);
% % xAxesVal = matObj.xAxesVal;
% % yAxesVal = matObj.yAxesVal;
% % axis([xAxesVal yAxesVal]);
% 
% % figure with patches
% 
% figure;
% colormap jet;
% % rivPtUnusedColor = [255,214,211]/255;
% % rivPtUnusedColor =[149, 144, 144]/255;
% % rivPtUnusedColor =[171, 168, 168]/255;
% rivPtUnusedColor =[180, 180, 180]/255;
% scatter([rivPt.xPos], [rivPt.yPos], 35, rivPtUnusedColor, 'filled');
% hold on;
% % scatter([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], 50, 'r', 'filled');
% delXText = 0.1;
% c = cellstr(num2str(rorOccupiedSet'));
% text([rivPt(rorOccupiedSet).xPos]+delXText, [rivPt(rorOccupiedSet).yPos], c, 'color', 'k','fontsize', 12);
% % colorbar;
% 
% % adding demand pts
% hold on;
% 
% for i = 1:length(rorOccupiedSet)
%     scatter(rivPt(rorOccupiedSet(i)).xPos, rivPt(rorOccupiedSet(i)).yPos, 100, colorVec{i}, 'filled','d');
%     hold on;
% end
% hold off;
% 
% matObj = matfile([pwd '/Data/demand/demandVars.mat']);
% xAxesVal = matObj.xAxesVal;
% yAxesVal = matObj.yAxesVal;
% axis([xAxesVal yAxesVal]);
% 
% % add patches
% xArr = [demPt.xPos] - demPt(1).xPos;
% nZXArr = xArr > 0;
% xDiff = min(xArr(nZXArr)); % closest neighbor in grid
% 
% yArr = [demPt.yPos] - demPt(1).yPos;
% nZYArr = yArr > 0;
% yDiff = min(yArr(nZYArr)); % closest neighbor in grid
% 
% hold on;
% for i = 1:length(D)
%     for j = 1:length(D{i})
%         demPtCoord = [demPt(D{i}(j)).xPos,demPt(D{i}(j)).yPos];
%         vert = repmat(demPtCoord, 4, 1) + [xDiff, yDiff; xDiff, -yDiff ; -xDiff, -yDiff; -xDiff, yDiff]/2;
%         patch('faces', [1,2,3,4], 'Vertices', vert, 'Facecolor', colorVec{i}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%     end
% end
% 
% xlabel('x (\times 100 km)');
% ylabel('y (\times 100 km)');
% titStr = sprintf('(numRivPt, numDemPt) = (%d, %d), \lambda = %4.3f', length(rivPt), length(demPt), lambda);
% title(titStr);
% grid;