% ====================================================================
% Annual TEOS-10 potential energy and annual conversion diagnostics.
% Use this file for updated annual potential energy. For event-stage conversion use RUN_eddy_fullpath_conversion.
% Reservoirs are in J/m^3; conversion rates are in W/m^3.
% ====================================================================

clc; clear; close all;

% Energy conversion diagnostics for the waters southwest of Taiwan Island.
% Every diagnostic is exported as an independent PNG. No subplot or title.
velocityFile = "D:\流浪C盘\CMEMS\2017的300米流速黑潮转向区.nc";
hydroFile = "D:\流浪C盘\CMEMS\2017的300米温盐黑潮转向区.nc";
outputRoot = "C:\Users\86182\Desktop\图片-去除subplot\Energy_conversion";
btDir = fullfile(outputRoot,"Barotropic_MKE_to_EKE");
bcDir = fullfile(outputRoot,"Baroclinic_MPE_to_EPE");
profileDir = fullfile(outputRoot,"Profiles");
energyRoot = fullfile(getenv('USERPROFILE'),'Desktop','图片-去除subplot','Energy');
mpeDir = fullfile(energyRoot,'MPE');
epeDir = fullfile(energyRoot,'EPE');
dirs = {outputRoot,btDir,bcDir,profileDir,energyRoot,mpeDir,epeDir};
for i = 1:numel(dirs)
    if ~exist(dirs{i},'dir'); mkdir(dirs{i}); end
end

lonRange = [115 123];
latRange = [19 23.5];
rho0 = 1025;
g = 9.81;
earthRadius = 6371000;
blockSize = 20;

lonAll = double(ncread(velocityFile,'longitude'));
latAll = double(ncread(velocityFile,'latitude'));
depth = double(ncread(velocityFile,'depth'));
time = ncread(velocityFile,'time');
iLon = find(lonAll >= lonRange(1) & lonAll <= lonRange(2));
iLat = find(latAll >= latRange(1) & latAll <= latRange(2));
lon = lonAll(iLon); lat = latAll(iLat);
nx = numel(lon); ny = numel(lat); nz = numel(depth); nt = numel(time);

% Pressure for TEOS-10 density calculation.
[lon2,lat2] = ndgrid(lon,lat);
pressure3 = zeros(nx,ny,nz,'single');
for k = 1:nz
    pressure3(:,:,k) = single(gsw_p_from_z(-depth(k).*ones(nx,ny),lat2));
end
[lon3,lat3] = ndgrid(lon,lat,depth);
% Outside the Baltic Sea, SA = SR(SP) + deltaSA(p,lon,lat).  Precomputing
% deltaSA is algebraically equivalent to a full conversion on every day.
deltaSA3 = gsw_SA_from_SP(35*ones(nx,ny,nz),double(pressure3),lon3,lat3) ...
    - gsw_SR_from_SP(35);

%% Pass 1: accumulate valid daily samples for the mean velocity and density fields.
sumU = zeros(nx,ny,nz); sumV = zeros(nx,ny,nz); sumR = zeros(nx,ny,nz);
sumSA = zeros(nx,ny,nz); sumCT = zeros(nx,ny,nz);
cntU = zeros(nx,ny,nz); cntV = zeros(nx,ny,nz); cntR = zeros(nx,ny,nz);
cntSA = zeros(nx,ny,nz); cntCT = zeros(nx,ny,nz);
for t0 = 1:blockSize:nt
    nb = min(blockSize,nt-t0+1);
    start = [iLon(1) iLat(1) 1 t0];
    count = [nx ny nz nb];
    u = clean_field(double(ncread(velocityFile,'uo',start,count)));
    v = clean_field(double(ncread(velocityFile,'vo',start,count)));
    temp = clean_field(double(ncread(hydroFile,'thetao',start,count)));
    sal = clean_field(double(ncread(hydroFile,'so',start,count)));
    pressure = repmat(double(pressure3),1,1,1,nb);
    SA = gsw_SR_from_SP(sal) + deltaSA3;
    CT = gsw_CT_from_pt(SA,temp);
    rho = gsw_rho(SA,CT,pressure);
    [sumU,cntU] = add_valid(sumU,cntU,u);
    [sumV,cntV] = add_valid(sumV,cntV,v);
    [sumR,cntR] = add_valid(sumR,cntR,rho);
    [sumSA,cntSA] = add_valid(sumSA,cntSA,SA);
    [sumCT,cntCT] = add_valid(sumCT,cntCT,CT);
    fprintf('Mean pass: %d/%d days\n',t0+nb-1,nt);
end
uMean = sumU./cntU; vMean = sumV./cntV; rhoMean = sumR./cntR;
uMean(cntU==0)=NaN; vMean(cntV==0)=NaN; rhoMean(cntR==0)=NaN;
SAmean = valid_mean(sumSA,cntSA); CTmean = valid_mean(sumCT,cntCT);
clear sumU sumV sumR sumSA sumCT cntU cntV cntR cntSA cntCT;

% Common annual reference column and TEOS-10 buoyancy frequency.
rhoRef = squeeze(mean(rhoMean,[1 2],'omitnan'));
SAref = squeeze(mean(SAmean,[1 2],'omitnan'));
CTref = squeeze(mean(CTmean,[1 2],'omitnan'));
pRef = squeeze(mean(double(pressure3),[1 2],'omitnan'));
[N2mid,pMid] = gsw_Nsquared(SAref(:),CTref(:),pRef(:),mean(lat));
N2 = interp1(pMid,N2mid,pRef,'linear','extrap');
N2(N2<=1e-7 | ~isfinite(N2)) = NaN;
N23 = reshape(N2,1,1,[]);
if sum(isfinite(N2)) < nz/2
    error('Too few valid N2 levels: %d/%d',sum(isfinite(N2)),nz);
end

%% Compute horizontal gradients of the mean fields per metre.
[dUdx,dUdy] = horizontal_gradient(uMean,lon,lat,earthRadius);
[dVdx,dVdy] = horizontal_gradient(vMean,lon,lat,earthRadius);
[dRdx,dRdy] = horizontal_gradient(rhoMean,lon,lat,earthRadius);

%% Pass 2: subtract the local temporal means and accumulate perturbation products.
sumUU=zeros(nx,ny,nz); sumVV=zeros(nx,ny,nz); sumUV=zeros(nx,ny,nz);
sumUR=zeros(nx,ny,nz); sumVR=zeros(nx,ny,nz);
sumRR=zeros(nx,ny,nz);
cntUU=zeros(nx,ny,nz); cntVV=zeros(nx,ny,nz); cntUV=zeros(nx,ny,nz);
cntUR=zeros(nx,ny,nz); cntVR=zeros(nx,ny,nz);
cntRR=zeros(nx,ny,nz);
for t0 = 1:blockSize:nt
    nb = min(blockSize,nt-t0+1);
    start = [iLon(1) iLat(1) 1 t0];
    count = [nx ny nz nb];
    u = clean_field(double(ncread(velocityFile,'uo',start,count)));
    v = clean_field(double(ncread(velocityFile,'vo',start,count)));
    temp = clean_field(double(ncread(hydroFile,'thetao',start,count)));
    sal = clean_field(double(ncread(hydroFile,'so',start,count)));
    pressure = repmat(double(pressure3),1,1,1,nb);
    SA = gsw_SR_from_SP(sal) + deltaSA3;
    CT = gsw_CT_from_pt(SA,temp);
    rho = gsw_rho(SA,CT,pressure);
    up = u-uMean; vp = v-vMean; rp = rho-rhoMean;
    [sumUU,cntUU]=add_valid(sumUU,cntUU,up.^2);
    [sumVV,cntVV]=add_valid(sumVV,cntVV,vp.^2);
    [sumUV,cntUV]=add_valid(sumUV,cntUV,up.*vp);
    [sumUR,cntUR]=add_valid(sumUR,cntUR,up.*rp);
    [sumVR,cntVR]=add_valid(sumVR,cntVR,vp.*rp);
    [sumRR,cntRR]=add_valid(sumRR,cntRR,rp.^2);
    fprintf('Covariance pass: %d/%d days\n',t0+nb-1,nt);
end
uu = valid_mean(sumUU,cntUU); vv = valid_mean(sumVV,cntVV);
uv = valid_mean(sumUV,cntUV); ur = valid_mean(sumUR,cntUR); vr = valid_mean(sumVR,cntVR);
rr = valid_mean(sumRR,cntRR);
clear sumUU sumVV sumUV sumUR sumVR sumRR cntUU cntVV cntUV cntUR cntVR cntRR;

%% BT: uu, uv and vv are window means of u_prime^2, u_prime*v_prime and v_prime^2.
C_BT = -rho0.*(uu.*dUdx + uv.*(dUdy+dVdx) + vv.*dVdy);

%% BC: ur and vr are mean horizontal perturbation density fluxes; N23 is the stratification factor on the grid.
C_BC = -(g^2./(rho0.*N23)).*(ur.*dRdx + vr.*dRdy);
MPE = g^2.*(rhoMean-reshape(rhoRef,1,1,[])).^2./(2*rho0.*N23);
EPE = g^2.*rr./(2*rho0.*N23);

%% Mask extreme values before exporting the fields.
C_BT = robust_mask(C_BT);
C_BC = robust_mask(C_BC);

targetDepths = [0 100 300];
depthLabels = {'0m','100m','300m'};
depthIndex = zeros(size(targetDepths));
for k = 1:numel(targetDepths)
    [~,depthIndex(k)] = min(abs(depth-targetDepths(k)));
end
btLim = common_limit(C_BT(:,:,depthIndex));
bcLim = common_limit(C_BC(:,:,depthIndex));
for k = 1:numel(depthIndex)
    label = sprintf('%dm',round(depth(depthIndex(k))));
    plot_map(C_BT(:,:,depthIndex(k)),lon,lat,btLim, ...
        fullfile(btDir,sprintf('CBT_%02d_%s.png',k,label)));
    plot_map(C_BC(:,:,depthIndex(k)),lon,lat,bcLim, ...
        fullfile(bcDir,sprintf('CBC_%02d_%s.png',k,label)));
    plot_energy_map(MPE(:,:,depthIndex(k)),lon,lat, ...
        fullfile(mpeDir,sprintf('MPE_%02d_%s.png',k,depthLabels{k})));
    plot_energy_map(EPE(:,:,depthIndex(k)),lon,lat, ...
        fullfile(epeDir,sprintf('EPE_%02d_%s.png',k,depthLabels{k})));
end

upper = find(depth<=150);
weights = layer_weights(depth(upper));
CBT150 = weighted_vertical_mean(C_BT(:,:,upper),weights);
CBC150 = weighted_vertical_mean(C_BC(:,:,upper),weights);
plot_map(CBT150,lon,lat,common_limit(CBT150),fullfile(btDir,'CBT_04_0-150m.png'));
plot_map(CBC150,lon,lat,common_limit(CBC150),fullfile(bcDir,'CBC_04_0-150m.png'));

btProfile = squeeze(mean(C_BT,[1 2],'omitnan'));
bcProfile = squeeze(mean(C_BC,[1 2],'omitnan'));
plot_profile(btProfile,depth,fullfile(profileDir,'CBT_area_mean_profile.png'));
plot_profile(bcProfile,depth,fullfile(profileDir,'CBC_area_mean_profile.png'));

summaryTable = table(depth(:),btProfile(:),bcProfile(:), ...
    'VariableNames',{'depth_m','CBT_W_m3','CBC_W_m3'});
writetable(summaryTable,fullfile(outputRoot,'area_mean_profiles.csv'));
save(fullfile(outputRoot,'energy_conversion_results.mat'), ...
    'lon','lat','depth','C_BT','C_BC','MPE','EPE','CBT150','CBC150','N2','-v7.3');
fprintf('N2 range: %.3g to %.3g s^-2; MPE/EPE medians: %.3g / %.3g J m^-3\n', ...
    min(N2,[],'omitnan'),max(N2,[],'omitnan'), ...
    median(MPE(:),'omitnan'),median(EPE(:),'omitnan'));

fprintf('Energy-conversion figures saved to %s\n',outputRoot);

function x = clean_field(x)
% Exclude missing or invalid input samples.
x(abs(x)>1e10)=NaN;
end

function [s,c] = add_valid(s,c,x)
% Accumulate valid sums and sample counts at each cell.
valid=isfinite(x);
x(~valid)=0;
s=s+sum(x,4);
c=c+sum(valid,4);
end

function m = valid_mean(s,c)
% Divide each sum by its local valid count.
m=s./c; m(c==0)=NaN;
end

function [dfdx,dfdy] = horizontal_gradient(f,lon,lat,R)
% Differentiate in the eastward and northward directions using physical distances.
[nx,ny,nz]=size(f); dfdx=nan(size(f)); dfdy=nan(size(f));
for k=1:nz
    for j=1:ny
        x=R*cosd(lat(j))*deg2rad(lon(:));
        dfdx(:,j,k)=gradient(f(:,j,k),x);
    end
    y=R*deg2rad(lat(:));
    for i=1:nx
        dfdy(i,:,k)=gradient(squeeze(f(i,:,k))',y)';
    end
end
end

function x = robust_mask(x)
% Apply the percentile-based extreme-value mask.
v=abs(x(isfinite(x)));
if isempty(v); return; end
q=prctile(v,99.9);
x(abs(x)>q)=NaN;
end

function lim = common_limit(x)
% Set the plotting limit from valid samples.
v=abs(x(isfinite(x)));
if isempty(v); lim=1; else; lim=prctile(v,98); end
if ~isfinite(lim) || lim<=0; lim=max(v); end
end

function w = layer_weights(z)
% Compute thickness weights for the selected levels.
z=z(:); e=[0;(z(1:end-1)+z(2:end))/2;z(end)+(z(end)-z(end-1))/2];
w=diff(e); w=w./sum(w);
end

function m = weighted_vertical_mean(x,w)
% Average across depth using layer thickness.
w=reshape(w,1,1,[]); valid=isfinite(x);
num=sum(x.*w,3,'omitnan'); den=sum(valid.*w,3);
m=num./den; m(den==0)=NaN;
end

function plot_map(field,lon,lat,lim,outputFile)
fig=figure('Visible','off','Color','w','Position',[100 80 820 680]);
m_proj('miller','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
m_pcolor(lon,lat,field'); shading interp; hold on;
caxis([-lim lim]); colormap(blue_white_red(256));
cb=colorbar; ylabel(cb,'W m^{-3}','FontSize',11);
m_coast('patch',[.82 .82 .82],'edgecolor',[.25 .25 .25]);
m_grid('fancy','tickdir','out');
m_text(120.7,22.7,'台湾岛','horizontalalignment','center', ...
    'verticalalignment','middle','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',10,'LineWidth',1);
exportgraphics(fig,outputFile,'Resolution',300); close(fig);
end

function plot_profile(x,depth,outputFile)
fig=figure('Visible','off','Color','w','Position',[100 80 620 700]);
plot(x,depth,'k-','LineWidth',2); hold on; xline(0,'--','Color',[.35 .35 .35]);
set(gca,'YDir','reverse','FontName','Times New Roman','FontSize',11,'LineWidth',1);
xlabel('Energy conversion (W m^{-3})'); ylabel('Depth (m)'); grid on; box on;
exportgraphics(fig,outputFile,'Resolution',300); close(fig);
end

function plot_energy_map(field,lon,lat,outputFile)
% Export one reservoir map.
fig=figure('Visible','off','Color','w','Position',[100 80 820 680]);
m_proj('miller','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
m_pcolor(lon,lat,field'); shading interp; hold on;
v=field(isfinite(field));
if isempty(v); lim=1; else; lim=prctile(v,98); end
caxis([0 max(lim,eps)]); colormap(parula(256));
cb=colorbar; ylabel(cb,'J m^{-3}','FontSize',11);
m_coast('patch',[.82 .82 .82],'edgecolor',[.25 .25 .25]);
m_grid('fancy','tickdir','out');
set(gca,'FontName','Times New Roman','FontSize',10,'LineWidth',1);
exportgraphics(fig,outputFile,'Resolution',300); close(fig);
end

function cmap = blue_white_red(n)
if nargin<1; n=256; end
x=[0 .5 1]; r=[0 1 .75]; g=[.25 1 0]; b=[.75 1 0];
xi=linspace(0,1,n); cmap=[interp1(x,r,xi)' interp1(x,g,xi)' interp1(x,b,xi)'];
end
