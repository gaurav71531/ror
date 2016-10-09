function getRiverData()
% riverData retrieval

format = 'html';
% siteLocation = {'siteID_Minnesota.txt', 'siteID_Iowa.txt', 'siteID_Illinois.txt', 'siteID_Missouri.txt', 'siteID_Wisconsin.txt'};
siteLocation = {'siteID_Nebraska.txt'};
stateID = {'ne'};
% stateID = {'mn', 'ia', 'il', 'mo', 'wi'};
header = [1,1,1,1,1,1];
% siteNum = '05378500';
% beginDate = '1991-01-01';
beginDate = '2008-01-01';
endDate = '2016-09-25';

for siteLocInd = 1:length(siteLocation)
    fName = sprintf('/Data/siteID/%s', siteLocation{siteLocInd});
    fid = fopen([pwd fName]);
    siteID = textscan(fid, '%s');
    fclose(fid);
    siteNum = siteID{1,1};
    for siteInd = header(siteLocInd):length(siteNum)
%         get data
        urlStr = sprintf('http://nwis.waterdata.usgs.gov/%s/nwis/uv?cb_00060=on&format=%s&site_no=%s&period=&begin_date=%s&end_date=%s', stateID{siteLocInd}, format, siteNum{siteInd}, beginDate, endDate);
%         urlStr = 'http://nwis.waterdata.usgs.gov/mn/nwis/uv?cb_00060=on&format=html&site_no=05387320&period=&begin_date=1991-01-01&end_date=2016-09-25';
        tableNum = 2;  %initial number for HTML table
        completeFlag = 0;
        numAttempt = 0;
        tic
        while(~completeFlag)
            try
                data = getTableFromWeb_mod(urlStr, tableNum);
            catch
%                 fprintf('flow table not available for siteID = %s from %s\n', siteNum{siteInd}, siteLocation{siteLocInd}(8:end-4));
%                 break;
            end
            if strcmp(data{1,1}, 'Date / Time')
                completeFlag = 1;
            end
            tableNum = tableNum + 1;
            numAttempt = numAttempt + 1;
            if numAttempt > 4
                fprintf('flow table read failed for siteID = %s from %s\n', siteNum{siteInd}, siteLocation{siteLocInd}(8:end-4));
                break;
            end
        end
        if completeFlag
            time = toc;
            fprintf('read succesful for siteID = %s from %s, time taken = %f\n', siteNum{siteInd}, siteLocation{siteLocInd}(8:end-4), time);
        else
            continue;
        end
%         save data
        %comment style used - #
        [nRow, ~] = size(data);
        data{1,1} = ['# ', data{1,1}];
        formatSpec = '%s %s\n';
        fName = sprintf('/Data/siteFlow/flow_siteID_%s.csv', siteNum{siteInd});
        fid = fopen([pwd fName], 'w');
        for row = 1:nRow
            fprintf(fid, formatSpec, data{row,:});
        end
        fclose(fid);
    end
end