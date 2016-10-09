function getPowerDemand()


csvLink = getCsvLinks();
saveRawData(csvLink);
saveElectricityUsageData();


function csvLink = getCsvLinks()

url = 'http://en.openei.org/datasets/files/961/pub/EPLUS_TMY2_RESIDENTIAL_BASE/';
html = urlread(url);
matches = regexp(html,'<a href=.*?/a>','match');
csvLink = cell(length(matches),1);
csvInd = 1;
% reqdStates = {'MN', 'IA', 'IL', 'WI', 'MO'};
reqdStates = {'AL','AK','AS','AZ','AR','CA','CO','CT','DE','DC','FL','GA','GU','HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MH','MA','MI','FM','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','MP','OH','OK','OR','PW','PA','PR','RI','SC','SD','TN','TX','UT','VT', 'VA', 'VI', 'WA', 'WV', 'WI', 'WY'};
for hyperInd = 1:length(matches)
    if strcmp(matches{hyperInd}(10:12), 'USA')
        for stateInd = 1:length(reqdStates)
            if strcmp(matches{hyperInd}(14:15), reqdStates{stateInd})
                reqdStr = regexp(matches{hyperInd}, '".*?"', 'match');
                csvLink{csvInd} = reqdStr{1}(2:end-1);
                csvInd = csvInd + 1;
                break;
            end
        end
    end
end
csvLink(csvInd:end) = [];


function saveRawData(csvLink)

timeOutVal = 10;
options = weboptions('Timeout',timeOutVal);
for dataInd = 1:length(csvLink)
    urlStr = sprintf('http://en.openei.org/datasets/files/961/pub/EPLUS_TMY2_RESIDENTIAL_BASE/%s', csvLink{dataInd}); 
    fileName = sprintf('/Data/demand/residentialUseRaw/dLoc_%d_%s.csv', dataInd, csvLink{dataInd}(5:6));
    if exist([pwd fileName], 'file') == 2,
        fprintf('skipping file dLoc_%d_%s.csv\n', dataInd, csvLink{dataInd}(5:6));
        continue;
    end
    nTry = 0;
    readSuccess = 0;
    while(nTry < 4 && ~readSuccess)
        try 
            outfilename = websave([pwd fileName],urlStr, options);
            readSuccess = 1;
        catch
            nTry = nTry + 1;
            options = weboptions('TimeOut', timeOutVal + 5*nTry);
        end
    end
    
    fprintf('data for dLoc_%d in %s saved\n', dataInd, csvLink{dataInd}(5:6));
end


function saveElectricityUsageData()

fileName = dir([pwd '/Data/demand/residentialUseRaw']);
for fileInd = 1:length(fileName)
    if length(fileName(fileInd).name) < 4, continue;end
    if ~strcmp(fileName(fileInd).name(end-2:end), 'csv'), continue;end
    fid = fopen([pwd '/Data/demand/residentialUseRaw/' fileName(fileInd).name]);
    data = textscan(fid, '%s %s %*[^\n]', 'delimiter', ',');
    dateAtime = data{1,1};
    powerDemand = data{1,2};
    fileNameCred = textscan(fileName(fileInd).name(1:end-4), '%s', 'delimiter', '_');
    fileStr = sprintf('dResiLoc_%s_%s.mat', fileNameCred{1}{2}, fileNameCred{1}{3});
    save([pwd '/Data/demand/elecPowDemand/' fileStr], 'dateAtime', 'powerDemand');
    fclose(fid);
end