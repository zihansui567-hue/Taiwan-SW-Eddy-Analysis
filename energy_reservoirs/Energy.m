% ====================================================================
% Annual kinetic energy and baseline potential energy fields (J/m^3).
% Use the kinetic energy outputs here and Energy_conversion.m for updated TEOS-10 potential energy.
% ====================================================================

clc; clear; close all;

% Each energy field/depth is exported as an independent PNG.
% No subplot and no title are used.
velocityFile = "D:\流浪C盘\CMEMS\2017的300米流速黑潮转向区.nc";
hydroFile = "D:\流浪C盘\CMEMS\2017的300米温盐黑潮转向区.nc";
outputRoot = fullfile(getenv('USERPROFILE'),'Desktop','图片-去除subplot','Energy');
fieldNames = {'MKE','EKE','MPE','EPE'};
for i = 1:numel(fieldNames)
    fieldDir = fullfile(outputRoot,fieldNames{i});
    if ~exist(fieldDir,'dir'); mkdir(fieldDir); end
end

lonRange = [115 123];
latRange = [19 23.5];
rho0 = 1025;
g = 9.81;

%% MKE and EKE
lon0 = double(ncread(velocityFile,'longitude'));
lat0 = double(ncread(velocityFile,'latitude'));
time = ncread(velocityFile,'time');
depth = double(ncread(velocityFile,'depth'));
iLon = find(lon0 >= lonRange(1) & lon0 <= lonRange(2));
iLat = find(lat0 >= latRange(1) & lat0 <= latRange(2));
lon = lon0(iLon); lat = lat0(iLat);
u = double(ncread(velocityFile,'uo',[iLon(1) iLat(1) 1 1], ...
    [numel(iLon) numel(iLat) numel(depth) numel(time)]));
v = double(ncread(velocityFile,'vo',[iLon(1) iLat(1) 1 1], ...
    [numel(iLon) numel(iLat) numel(depth) numel(time)]));

uMean = mean(u,4,'omitnan');
vMean = mean(v,4,'omitnan');
uPrime = u-uMean;
vPrime = v-vMean;

MKE4 = 0.5*rho0.*(uMean.^2+vMean.^2);
EKE4 = 0.5*rho0.*mean(uPrime.^2+vPrime.^2,4,'omitnan');

depthIndex = [1 min(22,numel(depth)) numel(depth)];
depthLabels = {'0m','100m','300m'};
for k = 1:3
    mkeMap = MKE4(:,:,depthIndex(k));
    ekeMap = EKE4(:,:,depthIndex(k));
    plot_energy_map(mkeMap,lon,lat,[0 300], ...
        fullfile(outputRoot,'MKE',sprintf('MKE_%02d_%s.png',k,depthLabels{k})));
    plot_energy_map(ekeMap,lon,lat,[0 80], ...
        fullfile(outputRoot,'EKE',sprintf('EKE_%02d_%s.png',k,depthLabels{k})));
end
clear u v uMean vMean uPrime vPrime MKE4 EKE4;

%% MPE and EPE
lon0 = double(ncread(hydroFile,'longitude'));
lat0 = double(ncread(hydroFile,'latitude'));
time = ncread(hydroFile,'time');
depth = double(ncread(hydroFile,'depth'));
iLon = find(lon0 >= lonRange(1) & lon0 <= lonRange(2));
iLat = find(lat0 >= latRange(1) & lat0 <= latRange(2));
lon = lon0(iLon); lat = lat0(iLat);
temp = double(ncread(hydroFile,'thetao',[iLon(1) iLat(1) 1 1], ...
    [numel(iLon) numel(iLat) numel(depth) numel(time)]));
sal = double(ncread(hydroFile,'so',[iLon(1) iLat(1) 1 1], ...
    [numel(iLon) numel(iLat) numel(depth) numel(time)]));

[~,Lat2] = ndgrid(lon,lat);
p3 = zeros(numel(lon),numel(lat),numel(depth));
for iz = 1:numel(depth)
    p3(:,:,iz) = gsw_p_from_z(-depth(iz).*ones(size(Lat2)),Lat2);
end
pressure = repmat(p3,1,1,1,numel(time));
[Lon3,Lat3] = ndgrid(lon,lat,depth);
deltaSA3 = gsw_SA_from_SP(35*ones(size(p3)),p3,Lon3,Lat3) ...
    - gsw_SR_from_SP(35);
SA = gsw_SR_from_SP(sal)+deltaSA3;
CT = gsw_CT_from_pt(SA,temp);
rho = gsw_rho(SA,CT,pressure);
clear temp sal pressure;

rhoRef = squeeze(mean(rho,[1 2 4],'omitnan'));
SAref = squeeze(mean(SA,[1 2 4],'omitnan'));
CTref = squeeze(mean(CT,[1 2 4],'omitnan'));
pRef = squeeze(mean(p3,[1 2],'omitnan'));
[N2mid,pMid] = gsw_Nsquared(SAref(:),CTref(:),pRef(:),mean(lat));
N2 = interp1(pMid,N2mid,pRef,'linear','extrap');
N2(N2 <= 1e-7 | ~isfinite(N2)) = NaN;
rhoRef4 = reshape(rhoRef,1,1,[],1);
N24 = reshape(N2,1,1,[],1);

rhoMean = mean(rho,4,'omitnan');
rhoPrime = rho-rhoMean;
MPE4 = g^2.*(rhoMean-rhoRef4).^2./(2*rho0.*N24);
EPE4 = g^2.*mean(rhoPrime.^2,4,'omitnan') ...
    ./(2*rho0.*N24);

depthIndex = [1 min(22,numel(depth)) numel(depth)];
for k = 1:3
    mpeMap = MPE4(:,:,depthIndex(k));
    epeMap = EPE4(:,:,depthIndex(k));
    plot_energy_map(mpeMap,lon,lat,[0 150], ...
        fullfile(outputRoot,'MPE',sprintf('MPE_%02d_%s.png',k,depthLabels{k})));
    plot_energy_map(epeMap,lon,lat,[0 400], ...
        fullfile(outputRoot,'EPE',sprintf('EPE_%02d_%s.png',k,depthLabels{k})));
end

fprintf('Energy single figures saved to %s\n',outputRoot);

function plot_energy_map(field,lon,lat,colorLimits,outputFile)
% Export one reservoir map.
fig = figure('Visible','off','Color','w','Position',[100 80 820 680]);
hold on;
m_proj('miller','lon',[115 123],'lat',[19 23.5]);
m_pcolor(lon,lat,field'); shading interp;
caxis(colorLimits); colormap(jet); colorbar;
m_coast('patch',[.8 .8 .8],'edgecolor',[.3 .3 .3]);
m_grid('fancy','tickdir','out');
m_text(120.7,22.7,'台湾岛','horizontalalignment','center', ...
    'verticalalignment','middle','FontSize',8,'FontWeight','bold');
exportgraphics(fig,outputFile,'Resolution',300);
close(fig);
end
