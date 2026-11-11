% ====================================================================
% Kuroshio South China Sea Index and counts by path class.
% Read daily ADT in metres and event indices from the spreadsheet.
% KSI has units m^2/s and includes the minus sign shown below.
% ====================================================================

clc; clear; close all;

% KSI = integral integral [-(g/f) * laplacian(eta)] dA
% Compute the Laplacian on the spherical grid before integrating over area.
adtFile   = "D:\流浪C盘\CMEMS_2014to2023_adt.nc";
eventFile = "D:\流浪C盘\SRDP\生成数量统计.xlsx";
outputDir = fullfile(getenv('USERPROFILE'), 'Desktop');
lonRange = [118.5 121];
latRange = [18.5 23];

lon0 = double(ncread(adtFile, 'longitude'));
lat0 = double(ncread(adtFile, 'latitude'));
timeRaw = double(ncread(adtFile, 'time'));
timeDT = datetime(1970,1,1) + seconds(timeRaw);

iLon = find(lon0 >= lonRange(1) & lon0 <= lonRange(2));
iLat = find(lat0 >= latRange(1) & lat0 <= latRange(2));
lon = lon0(iLon);
lat = lat0(iLat);
eta = double(ncread(adtFile, 'adt', ...
    [iLon(1) iLat(1) 1], [numel(iLon) numel(iLat) numel(timeRaw)]));

KSI = compute_ksi(eta, lon, lat);
mu = mean(KSI, 'omitnan');
sigma = std(KSI, 0, 'omitnan');
fprintf('KSI mean = %.6e m^2 s^-1\n', mu);
fprintf('KSI std  = %.6e m^2 s^-1\n', sigma);

% Read event indices. Looping exceeds mean + SD; leaping falls below mean - SD.
eventData = readmatrix(eventFile);
atime = sanitize_indices(eventData(:,8), numel(KSI));
ctime = sanitize_indices(eventData(:,9), numel(KSI));

fig1 = figure('Visible','off','Color','w','Position',[80 80 1450 420]);
hKSI = plot(timeDT, KSI, 'k-', 'LineWidth', 0.9); hold on; grid on; box on;
hMu = yline(mu, 'k-', '\mu', 'LineWidth', 0.8);
hUp = yline(mu + sigma, 'k--', '\mu+\sigma', 'LineWidth', 0.8);
hDn = yline(mu - sigma, 'k--', '\mu-\sigma', 'LineWidth', 0.8);
hACE = scatter(timeDT(atime), KSI(atime), 18, 'o', ...
    'MarkerEdgeColor',[0.85 0 0], 'MarkerFaceColor','none', 'LineWidth',0.8);
hCE = scatter(timeDT(ctime), KSI(ctime), 18, 'o', ...
    'MarkerEdgeColor',[0 0.25 0.9], 'MarkerFaceColor','none', 'LineWidth',0.8);
xlabel('Time'); ylabel('KSI (m^2 s^{-1})');
title('Kuroshio South China Sea Index');
legend([hKSI hMu hUp hDn hACE hCE], ...
    {'KSI','\mu','\mu+\sigma','\mu-\sigma','ACE','CE'}, ...
    'Location','southoutside','Orientation','horizontal');
exportgraphics(fig1, fullfile(outputDir,'KSI_time_series.png'), 'Resolution',300);

% Count events in each threshold-defined path class.
Acounts = classify_counts(KSI(atime), mu, sigma);
Ccounts = classify_counts(KSI(ctime), mu, sigma);
fprintf('ACE [Looping, Leaking, Leaping] = [%d, %d, %d]\n', Acounts);
fprintf('CE  [Looping, Leaking, Leaping] = [%d, %d, %d]\n', Ccounts);

fig2 = figure('Visible','off','Color','w','Position',[120 120 760 520]);
cats = categorical({'Looping','Leaking','Leaping'}, ...
    {'Looping','Leaking','Leaping'});
vals = [Acounts(:), Ccounts(:)];
b = bar(cats, vals); grid on; box on;
b(1).FaceColor = [0.85 0.1 0.1];
b(2).FaceColor = [0.1 0.3 0.9];
ylabel('Number of eddies');
legend({'ACE','CE'},'Location','northoutside','Orientation','horizontal');
for k = 1:numel(b)
    text(b(k).XEndPoints, b(k).YEndPoints, string(b(k).YData), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
exportgraphics(fig2, fullfile(outputDir,'KSI_path_counts.png'), 'Resolution',300);

% Compare the 2017 empirical KSI density with a fitted Gaussian.
is2017 = year(timeDT) == 2017 & isfinite(KSI);
ksi2017 = KSI(is2017);
if isempty(ksi2017)
    error('KSI_analysis:No2017Data', '2017 年没有可用 KSI 数据。');
end
mu2017 = mean(ksi2017);
sigma2017 = std(ksi2017);
xFit = linspace(min(ksi2017), max(ksi2017), 300);
pdfFit = exp(-0.5*((xFit-mu2017)/sigma2017).^2) ./ (sigma2017*sqrt(2*pi));

fig3 = figure('Visible','off','Color','w','Position',[120 120 900 480]);
histogram(ksi2017,'Normalization','pdf','BinMethod','fd', ...
    'FaceColor',[0.3 0.6 0.8],'EdgeColor','k','FaceAlpha',0.8); hold on;
plot(xFit,pdfFit,'k-','LineWidth',1.4); grid on; box on;
xlabel('KSI (m^2 s^{-1})'); ylabel('Probability density');
title('2017 KSI distribution');
legend({'Histogram','Normal fit'},'Location','best');
exportgraphics(fig3, fullfile(outputDir,'KSI_2017_distribution.png'), 'Resolution',300);

% Summarize the full record and test normality.
ksiValid = KSI(isfinite(KSI));
n = numel(ksiValid);
muAll = mean(ksiValid);
sigmaAll = std(ksiValid);
alpha = 0.05;
zAlpha = -sqrt(2)*erfcinv(2*(1-alpha/2));
muCI = muAll + [-1 1] * zAlpha * sigmaAll/sqrt(n);
dataCI = prctile(ksiValid,[2.5 97.5]);
if exist('kstest','file') == 2 && sigmaAll > 0
    [hKS,pKS,ksStat] = kstest((ksiValid-muAll)/sigmaAll);
else
    hKS = NaN; pKS = NaN; ksStat = NaN;
end
fprintf('95%% mean CI = [%.6e, %.6e] m^2 s^-1\n',muCI(1),muCI(2));
fprintf('95%% empirical range = [%.6e, %.6e] m^2 s^-1\n',dataCI(1),dataCI(2));
fprintf('KS test: h=%g, p=%.4g, statistic=%.4g\n',hKS,pKS,ksStat);

xAll = linspace(min(ksiValid),max(ksiValid),300);
pdfAll = exp(-0.5*((xAll-muAll)/sigmaAll).^2) ./ (sigmaAll*sqrt(2*pi));
fig4 = figure('Visible','off','Color','w','Position',[120 120 900 500]);
hHist = histogram(ksiValid,30,'Normalization','pdf', ...
    'FaceColor',[0.3 0.6 0.8],'EdgeColor','k','FaceAlpha',0.8); hold on;
hFit = plot(xAll,pdfAll,'k-','LineWidth',1.4);
yMax = max(ylim);
hRange = plot(dataCI,[0.04 0.04]*yMax,'r-','LineWidth',2);
plot(dataCI,[0.04 0.04]*yMax,'r|','MarkerSize',15,'LineWidth',1.5);
grid on; box on;
xlabel('KSI (m^2 s^{-1})'); ylabel('Probability density');
title('KSI distribution and 95% empirical range');
legend([hHist hFit hRange],{'Histogram','Normal fit','95% data range'},'Location','best');
if isfinite(pKS)
    text(0.02,0.97,sprintf('KS p = %.3g',pKS),'Units','normalized', ...
        'VerticalAlignment','top','BackgroundColor','w');
end
exportgraphics(fig4, fullfile(outputDir,'KSI_confidence_interval.png'), 'Resolution',300);

fprintf('KSI figures saved to %s\n', outputDir);

function idx = sanitize_indices(raw, nMax)
% Retain valid event indices within the available time series.
idx = unique(round(raw(isfinite(raw))),'stable');
idx = idx(idx >= 1 & idx <= nMax);
end

function counts = classify_counts(values, mu, sigma)
% Classify event-day KSI relative to the supplied mean and standard deviation.
values = values(isfinite(values));
counts = [sum(values > mu+sigma), ...
          sum(values >= mu-sigma & values <= mu+sigma), ...
          sum(values < mu-sigma)];
end

function KSI = compute_ksi(eta, lon, lat)
% Integrate the signed vorticity estimate using grid-cell areas.
g = 9.81;
omega = 7.292115e-5;
Re = 6371000;
lon = double(lon(:));
lat = double(lat(:));
phi = deg2rad(lat);
lambda = deg2rad(lon);
dLambda = mean(diff(lambda));
dPhi = mean(diff(phi));
nx = numel(lon); nt = size(eta,3);

latEdges = [phi(1)-0.5*dPhi; 0.5*(phi(1:end-1)+phi(2:end)); phi(end)+0.5*dPhi];
cellAreaLat = Re^2 * abs(dLambda) .* abs(diff(sin(latEdges)));
dA = repmat(cellAreaLat.', nx, 1);
f = 2*omega*sin(phi).';

phiC = phi(2:end-1).';
cosC = cos(phiC);
cosNorth = cos(0.5*(phi(2:end-1)+phi(3:end))).';
cosSouth = cos(0.5*(phi(1:end-2)+phi(2:end-1))).';
KSI = nan(nt,1);

for t = 1:nt
    e = double(eta(:,:,t));
    eC = e(2:end-1,2:end-1);
    lonTerm = (e(3:end,2:end-1)-2*eC+e(1:end-2,2:end-1)) ...
        ./ (dLambda^2 .* cosC.^2);
    fluxNorth = cosNorth .* (e(2:end-1,3:end)-eC) ./ dPhi;
    fluxSouth = cosSouth .* (eC-e(2:end-1,1:end-2)) ./ dPhi;
    latTerm = (fluxNorth-fluxSouth) ./ (dPhi .* cosC);
    lap = (lonTerm + latTerm) ./ Re^2;
    integrand = -(g ./ f(2:end-1)) .* lap;
    weights = dA(2:end-1,2:end-1);
    valid = isfinite(integrand) & isfinite(weights);
    if any(valid,'all')
        KSI(t) = sum(integrand(valid).*weights(valid));
    end
end
end
