clear
close all;
% % % % % 
% month = march

dBegin = [3,0,0];

fileName = dir([pwd '/Data/demand/elecPowDemand/']);
reqdStates = {'MN', 'IA', 'IL', 'WI', 'MO'};
% mu = zeros(2,length(fileName));
% sigSq = zeros(2,length(fileName));
% alpha = zeros(1,length(fileName));
matObj = matfile([pwd '/Data/demand/demandVars.mat']);
mu = matObj.mu;
sigSq = matObj.sigSq;
alpha = matObj.alpha;
paramInd = 1;

for fileInd = 1:length(fileName)
% for fileInd = 5:10
    if length(fileName(fileInd).name) < 4, continue;end
    if ~strcmp(fileName(fileInd).name(end-2:end), 'mat'), continue;end
    fileAcceptFlag = 0;
    for stateInd = 1:length(reqdStates)
        nameUse = textscan(fileName(fileInd).name(1:end-4), '%s', 'delimiter', '_');
        if strcmp(nameUse{1}{3}, reqdStates{stateInd}),
            fileAcceptFlag = 1;
            break;
        end
    end
    if ~fileAcceptFlag, continue;end
    matObj = matfile([pwd '/Data/demand/elecPowDemand/' fileName(fileInd).name]);
    dateAtime = matObj.dateAtime;
    powerDemand = matObj.powerDemand;
    
    powVal = zeros(length(dateAtime),1);
    powInd = 1;
    for dataInd = 1:length(dateAtime)
        if str2double(dateAtime{dataInd}(1:2)) ==  dBegin(1),
            powVal(powInd) = str2double(powerDemand{dataInd});
            powInd = powInd + 1;
        end
    end
    powVal(powInd:end) = [];
%     [f,xi] = ksdensity(powVal);
    figure;
%     
%     plot(xi, f);grid;
%     xlabel('hourly power demand (kWh)');
%     ylabel('pdf');
%     figure;
%     
% %     hold on;
%     [counts, centers] = hist(powVal, 100);
%     
%     stem(centers, counts/(powInd-1));
%     
%     mu = mean(powVal);
%     sigSq = var(powVal);
%     
    x = linspace(min(powVal), max(powVal), 100);
    deltaX = mean(diff(x));
    GMModel = fitgmdist(powVal,2);
    y = pdf(GMModel,x');
    
    muUse_1 = mu.var(1) * randn() + mu.mean(1);
    muUse_2 = mu.var(2) * randn() + mu.mean(2);
%     sigSqUse_1 = sigSq.new(1) * sigSq.sigSq(1) / chi2rnd(sigSq.new(1));
%     sigSqUse_2 = sigSq.new(2) * sigSq.sigSq(2) / chi2rnd(sigSq.new(2));
    ygg = alpha * 1/sqrt(2*pi*sigSq(1))*exp(-(x-muUse_1).^2/2/sigSq(1)) + (1-alpha) * 1/sqrt(2*pi*sigSq(2))*exp(-(x-muUse_2).^2/2/sigSq(2));
    plot(x, y*deltaX);grid;
    hold on;
    plot(x, ygg*deltaX);grid;
end
% 
% mu(:,paramInd:end) = [];
% sigSq(:,paramInd:end) = [];
% alpha(paramInd:end) = [];
% 
% figure;
% plot(1:paramInd-1,mu(1,:), 1:paramInd-1,mu(2,:));grid;
% 
% figure;
% plot(1:paramInd-1, sigSq(1,:), 1:paramInd-1, sigSq(2,:));grid;
% 
% figure;
% plot(1:paramInd-1, alpha);grid;

% var = [0.022, 0.044];
% mu = [0.92, 1.7];
% alpha = 0.25;

% mu.mean = [0.92, 1.7];
% mu.var = [0.1, 0.1];
% sigSq.new = [0.022, 0.044];
% sigSq.sigSq = [0.1, 0.1];