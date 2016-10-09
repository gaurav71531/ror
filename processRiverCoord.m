clear;
% close all;
dataFile = {'latLon_Minnesota.mat', 'latLon_Iowa.mat', 'latLon_Illinois.mat', 'latLon_Missouri.mat', 'latLon_Wisconsin.mat'};
% dataFile = {'latLon_Minnesota.mat'};
figure;
colorVec = {'b', 'r', 'g', 'm', 'c', 'k'};
legendStr = cell(1,length(dataFile));
for i = 1:length(dataFile)
    load([pwd '/Data/siteLatLon/' dataFile{i}]);
    dczone = utmzone(mean(latLon(:,1),'omitnan'),mean(latLon(:,2),'omitnan'));
    utmstruct = defaultm('utm');
    utmstruct.zone = dczone;
    utmstruct.zone = '15T'; % plot multi-zone data in single zone for map consistency
    utmstruct.geoid = wgs84Ellipsoid('meters');
    
    utmstruct = defaultm(utmstruct);
    [x,y] = mfwdtran(utmstruct,latLon(:,1),latLon(:,2));
%     [x,y,~] = deg2utm(latLon(:,1),latLon(:,2));

    hold on;
    scattStr = sprintf('s%d = scatter(x, y, colorVec{i});', i);
    eval(scattStr);
%      scatter bug 1283854  %
    scattStr = sprintf('s%d.MarkerEdgeColor = s%d.CData;', i,i);
    eval(scattStr);
    
    a = [1:length(x)]';
    c = cellstr(num2str(a));

    text(x, y, c,'color', colorVec{i},'fontsize', 15);
    legStr = dataFile{i}(8:end-4);
    legendStr(i) = {legStr};
end
grid;
legend(legendStr);