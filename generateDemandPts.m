clear

matObj = matfile([pwd '/Data/demand/demandVarsMar.mat']);

numDemandPts = matObj.numDemandPts;
xMin = matObj.xMin;
xMax = matObj.xMax;
yMin = matObj.yMin;
yMax = matObj.yMax;

[x_d, y_d] = meshgrid(linspace(xMin, xMax, numDemandPts(1)), linspace(yMin, yMax, numDemandPts(2)));
x_d = x_d';
y_d = y_d';
x_d = x_d(:);
y_d = y_d(:);

% matObj = matfile([pwd '/Data/demand/demandVarsMar.mat']);
% matObj = matfile([pwd '/Data/demand/demandVarsJun.mat']);
matObj = matfile([pwd '/Data/demand/demandVarsNov.mat']);
mu = matObj.mu;
sigSq = matObj.sigSq;
alpha = matObj.alpha;

muUse_1 = mu.var(1) * randn(prod(numDemandPts),1) + mu.mean(1);
muUse_2 = mu.var(2) * randn(prod(numDemandPts),1) + mu.mean(2);

alpha = 0.05 * randn(prod(numDemandPts),1) + alpha;
demPt.x = x_d;
demPt.y = y_d;
demPt.alpha = alpha;
% demPt.alpha = ones(1,length(x_d))*alpha;
demPt.mu1 = muUse_1;
demPt.mu2 = muUse_2;
demPt.sigSq1 = ones(length(x_d),1)*sigSq(1);
demPt.sigSq2 = ones(length(x_d),1)*sigSq(2);

% for i = 1:prod(numDemandPts)
%     demPt(i).x = x_d(i);
%     demPt(i).y = y_d(i);
%     demPt(i).alpha = alpha;
%     demPt(i).mu1 = muUse_1(i);
%     demPt(i).mu2 = muUse_2(i);
%     demPt(i).sigSq1 = sigSq(1);
%     demPt(i).sigSq2 = sigSq(2);
% %     demPt(i).dPdf = @(x) alpha * 1/sqrt(2*pi*sigSq(1))*exp(-(x-muUse_1(i)).^2/2/sigSq(1)) + (1-alpha) * 1/sqrt(2*pi*sigSq(2))*exp(-(x-muUse_2(i)).^2/2/sigSq(2));
% end

% save([pwd '/Data/demand/demandProfileMar.mat'], 'demPt');
% save([pwd '/Data/demand/demandProfileJun.mat'], 'demPt');
save([pwd '/Data/demand/demandProfileNov.mat'], 'demPt');


