function getRiverCoord(varargin)

%coordinates
% % % % % % % % get side IDs
% e.g. getRiverCoord([pwd '/Data/siteID/siteID_Nebraska.txt'])
if nargin < 1,
    siteNum = {'05355092'}; 
else
    fid = fopen(varargin{1});
    siteID = textscan(fid, '%s');
    fclose(fid);
    siteNum = siteID{1,1};
end

% [latLon, meanMedian] = getLatLonMeanMedianMatrix(siteNum);
saveLatLon(siteNum);


function saveLatLon(siteNum)

for siteInd = 1:length(siteNum)
    
    urlStr = sprintf('http://waterdata.usgs.gov/nwis/inventory/?site_no=%s&agency_cd=USGS', siteNum{siteInd});
    data = webread(urlStr);

    k = strfind(data, 'Latitude');
    k1 = strfind(data(k:k+30), '"');
    latStr = data(k+8+2:k+k1-1);
%     latStr = strrep(latStr, '&#176;', 'd');
    latStr = textscan(latStr, '%s', 'Delimiter', {'&#176;','''', '"'});
    latVal = str2double({latStr{1}{1}, latStr{1}{2}, latStr{1}{3}});
    latVal = dms2degrees(latVal);

    k = strfind(data, 'Longitude');
    k1 = strfind(data(k:k+30), '"');
    lonStr = data(k+9+2:k+k1-1);
%     lonStr = strrep(lonStr, '&#176;', 'd');
    lonStr = textscan(lonStr, '%s', 'Delimiter', {'&#176;','''', '"'});
    lonVal = str2double({lonStr{1}{1}, lonStr{1}{2}, lonStr{1}{3}});
    lonVal(1) = -lonVal(1); % to take care of the W
    lonVal = dms2degrees(lonVal);
    fileStr = sprintf('/Data/siteLatLon/siteID_Loc_%s.txt', siteNum{siteInd});
    fid = fopen([pwd fileStr], 'w');
    fprintf(fid, '%f %f\n', latVal, lonVal);
    fclose(fid);
    fprintf('read data for siteID = %s\n', siteNum{siteInd});
end


function getLatLonMeanMedianMatrix(siteNum)


latLon = zeros(length(siteNum), 2);
meanMedian = zeros(length(siteNum), 2);
for siteInd = 1:length(siteNum)
    
    urlStr = sprintf('http://waterdata.usgs.gov/nwis/inventory/?site_no=%s&agency_cd=USGS', siteNum{siteInd});
    data = webread(urlStr);

    k = strfind(data, 'Latitude');
    k1 = strfind(data(k:k+30), '"');
    latStr = data(k+8+2:k+k1-1);
%     latStr = strrep(latStr, '&#176;', 'd');
    latStr = textscan(latStr, '%s', 'Delimiter', {'&#176;','''', '"'});
    latVal = str2double({latStr{1}{1}, latStr{1}{2}, latStr{1}{3}});
    latLon(siteInd,1) = dms2degrees(latVal);

    k = strfind(data, 'Longitude');
    k1 = strfind(data(k:k+30), '"');
    lonStr = data(k+9+2:k+k1-1);
%     lonStr = strrep(lonStr, '&#176;', 'd');
    lonStr = textscan(lonStr, '%s', 'Delimiter', {'&#176;','''', '"'});
    lonVal = str2double({lonStr{1}{1}, lonStr{1}{2}, lonStr{1}{3}});
    lonVal(1) = -lonVal(1); % to take care of the W
    latLon(siteInd,2) = dms2degrees(lonVal);
    
    urlStr = sprintf('http://waterdata.usgs.gov/mn/nwis/uv/?site_no=%s&PARAmeter_cd=00065,00060', siteNum{siteInd});
    tableNum = 3;  %initial number for HTML table
    completeFlag = [0,0];
    numAttempt = 0;
    while(~prod(completeFlag))
        try
            data = getTableFromWeb_mod(urlStr, tableNum);
        catch
            fprintf('Mean Median table not available for siteID = %s\n', siteNum{siteInd});
            break;
        end
        for i = 1:size(data,2)
            if strcmp(data{1,i}, 'Median'),
                meanMedian(siteInd, 2) = str2double(data{2,i});
                completeFlag(1) = 1;
            end
            if strcmp(data{1,i}, 'Mean'),
                meanMedian(siteInd, 1) = str2double(data{2,i});
                completeFlag(2) = 1;
            end
        end
        tableNum = tableNum + 1;
        numAttempt = numAttempt + 1;
        if numAttempt > 4,
            fprintf('Mean Median table read failed for siteID = %s\n', siteNum{siteInd});
            break;
        end
    end
    fprintf('read data for siteID = %s\n', siteNum{siteInd});
    assignin('base', 'latLon', latLon);
    assignin('base', 'meanMedian', meanMedian);
end

% processRiverCoord(latLon, meanMedian);


% function processRiverCoord(latLon, meanMedian)
% 
% dczone = utmzone(mean(latLon(:,1),'omitnan'),mean(latLon(:,2),'omitnan'));
% utmstruct = defaultm('utm');
% utmstruct.zone = dczone;
% utmstruct.geoid = wgs84Ellipsoid;
% utmstruct = defaultm(utmstruct);
% [x,y] = mfwdtran(utmstruct,latLon(:,1),latLon(:,2));
% 
% figure;
% 
% scatter(x, y);grid;
% 
% dx = 0.1 * mean(x);
% dy = 0.1 * mean(y);
% 
% c = cellstr(1:length(x));
% 
% text(x+dx, y+dy, c);



