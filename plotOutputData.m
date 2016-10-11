function plotOutputData(varargin)


colorVec = {[1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1], [0.7 0.1 0.9], [0 0 0], [207 206 102]/255, [82 172 202]/255};
% colorVec = {'b', 'r', 'g', 'm', 'c', 'k', 'y', 'w'};

% currentData is original without std
% currentDataNew is with some modifications in pdf mean and std
% currentDataNew1.mat is with Nebraska
% currentDataNewTest.mat is with fft values

% matObj = matfile('currentDataNew1.mat');
matObj = matfile('currentDataJun.mat');
rivPt = matObj.rivPt;

if nargin < 2
    lambda = 0.01;
    inpFileStr = sprintf('output_lambda_%4.3f_numDemPt_576_Jun.mat', lambda);
    matObj = matfile([pwd '/OutputData/' inpFileStr]);
    D = matObj.D;
    rorOccupiedSet = matObj.rorOccupiedSet;
    demPt = matObj.demPt;
    rivPt = matObj.rivPt;
%     lambda = matObj.lambda;
else
    D = varargin{1};
    rorOccupiedSet = varargin{2};
    demPt = varargin{3};
    lambda = varargin{4};
end
D = D(1:length(rorOccupiedSet));
    

% reorder based on cardinality
cardSet = cellfun(@(x) length(x), D);
[~, ind] = sort(cardSet, 'descend');
D = D(ind);
rorOccupiedSet = rorOccupiedSet(ind);

% figure with patches

figure;
colormap jet;
% rivPtUnusedColor = [255,214,211]/255;
% rivPtUnusedColor =[149, 144, 144]/255;
% rivPtUnusedColor =[171, 168, 168]/255;
rivPtUnusedColor =[180, 180, 180]/255;
scatter([rivPt.xPos], [rivPt.yPos], 35, rivPtUnusedColor, 'filled');
hold on;
% scatter([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], 50, 'r', 'filled');

matObj = matfile([pwd '/Data/demand/demandVars.mat']);
xAxesVal = matObj.xAxesVal;
yAxesVal = matObj.yAxesVal;
axis([xAxesVal yAxesVal]);

% add patches
xMin = min([demPt.xPos]);
xArr = [demPt.xPos] - xMin;
nZXArr = xArr > 0;
xDiff = min(xArr(nZXArr)); % closest neighbor in grid

yMin = min([demPt.yPos]);
yArr = [demPt.yPos] - yMin;
nZYArr = yArr > 0;
yDiff = min(yArr(nZYArr)); % closest neighbor in grid

xLen = floor((max([demPt.xPos]) - xMin) / xDiff) + 1;
yLen = floor((max([demPt.yPos]) - yMin) / yDiff) + 1;

demPtRORMatrix = zeros(xLen, yLen);
for rorInd = 1:length(D)
    for dInd = D{rorInd}
        demPtRORMatrix(floor((demPt(dInd).xPos-xMin)/xDiff)+1,floor((demPt(dInd).yPos-yMin)/yDiff)+1) = rorOccupiedSet(rorInd);
    end
end
alphaValFace = 0.15;
alphaValEdge = 0;

hold on;
for i = 1:length(D)
    for j = 1:length(D{i})
        demPtCoord = [demPt(D{i}(j)).xPos,demPt(D{i}(j)).yPos];
        vert = repmat(demPtCoord, 4, 1) + [xDiff, yDiff; xDiff, -yDiff ; -xDiff, -yDiff; -xDiff, yDiff]/2;
        
        xInd = floor((demPt(D{i}(j)).xPos-xMin)/xDiff)+1;
        yInd = floor((demPt(D{i}(j)).yPos-yMin)/yDiff)+1;
        faceColor = zeros(4,3);
        edgeAlpha = ones(4,1);
        
%         facecolor 1
        if xInd == xLen
            faceColor(1,:) = colorVec{i};
            edgeAlpha(1) = alphaValEdge;
        elseif demPtRORMatrix(xInd,yInd) == demPtRORMatrix(xInd+1,yInd)
            faceColor(1,:) = colorVec{i};
            edgeAlpha(1) = alphaValEdge;
        end
        
%         facecolor 2
        if yInd == 1
            faceColor(2,:) = colorVec{i};
            edgeAlpha(2) = alphaValEdge;
        elseif demPtRORMatrix(xInd,yInd) == demPtRORMatrix(xInd,yInd-1)
            faceColor(2,:) = colorVec{i};
            edgeAlpha(2) = alphaValEdge;
        end
        
%         facecolor 3
        if xInd == 1
            faceColor(3,:) = colorVec{i};
            edgeAlpha(3) = alphaValEdge;
        elseif demPtRORMatrix(xInd,yInd) == demPtRORMatrix(xInd-1,yInd)
            faceColor(3,:) = colorVec{i};
            edgeAlpha(3) = alphaValEdge;
        end
        
%         facecolor 4
        if yInd == yLen
            faceColor(4,:) = colorVec{i};
            edgeAlpha(4) = alphaValEdge;
        elseif demPtRORMatrix(xInd,yInd) == demPtRORMatrix(xInd,yInd+1)
            faceColor(4,:) = colorVec{i};
            edgeAlpha(4) = alphaValEdge;
        end
%         edgeAlpha = ones(4,1);
        
%         patch('faces', [1,2,3,4], 'Vertices', vert, 'Facecolor', colorVec{i}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        patch('faces', [1,2,3,4], 'Vertices', vert, 'Facecolor', colorVec{i},...
            'FaceAlpha', alphaValFace, 'EdgeColor', 'flat',...
            'FaceVertexCData', faceColor, 'EdgeAlpha', 'flat',...
            'FaceVertexAlphaData', edgeAlpha,...
            'AlphaDataMapping', 'none');
    end
end

delXText = 0.1;
% c = cellstr(num2str(rorOccupiedSet'));
c = cellstr(num2str([1:length(rorOccupiedSet)]'));
text([rivPt(rorOccupiedSet).xPos]+delXText, [rivPt(rorOccupiedSet).yPos], c, 'color', 'k','fontsize', 12, 'FontWeight', 'bold');
% colorbar;

% adding demand pts
hold on;

for i = 1:length(rorOccupiedSet)
    scatter(rivPt(rorOccupiedSet(i)).xPos, rivPt(rorOccupiedSet(i)).yPos, 150, colorVec{i}, 'filled','d');
    hold on;
end

xlabel('$x (\times 100\hspace{2pt}km)$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('$y (\times 100\hspace{2pt}km)$', 'interpreter', 'latex', 'FontSize', 14);
titStr = sprintf('$\\vert\\mathcal{R}\\vert$ = %d, $\\vert\\mathcal{D}\\vert$ = %d, $\\lambda$ = %4.3f', length(rivPt), length(demPt), lambda);
title(titStr, 'interpreter', 'latex', 'FontSize', 14);
% set(gca, 'TickLength', [0,0]);
box on;
% grid;

% figure with dots only

% figure;
% colormap jet;
% scatter([rivPt.xPos], [rivPt.yPos], 35, [255,214,211]/255, 'filled');
% hold on;
% % scatter([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], 50, 'r', 'filled');
% c = cellstr(num2str(rorOccupiedSet'));
% text([rivPt(rorOccupiedSet).xPos], [rivPt(rorOccupiedSet).yPos], c, 'color', 'k','fontsize', 12);
% % colorbar;
% xlabel('x (\times 100 km)');
% ylabel('y (\times 100 km)');
% grid;
% 
% % adding demand pts
% hold on;
% 
% for i = 1:length(rorOccupiedSet)
%     scatter(rivPt(rorOccupiedSet(i)).xPos, rivPt(rorOccupiedSet(i)).yPos, 70, colorVec{i}, 'filled');
%     hold on;
%     scatter([demPt(D{i}).xPos], [demPt(D{i}).yPos], 10, colorVec{i});
%     hold on;
% end
% hold off;
% 
% matObj = matfile([pwd '/Data/demand/demandVars.mat']);
% xAxesVal = matObj.xAxesVal;
% yAxesVal = matObj.yAxesVal;
% axis([xAxesVal yAxesVal]);