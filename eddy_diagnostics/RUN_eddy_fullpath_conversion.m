function RUN_eddy_fullpath_conversion
% ====================================================================
% Barotropic and baroclinic conversion across eight eddy stages.
% Each target represents the 30-day window from day -15 through day +14.
% Positive BT transfers MKE to EKE; positive BC transfers mean to eddy available potential energy (W/m^3).
% ====================================================================

% Stage diagnostics. 30-day local Reynolds decomposition, NOT daily instantaneous fields.
root='C:/Users/86182/Desktop/图片-去除subplot/Eddy_fullpath_20261106';
if ~isfolder(root), mkdir(root); end
diary(fullfile(root,'run_log.txt')); done=onCleanup(@()diary('off'));
addpath(genpath('C:/Users/86182/Desktop/gsw_matlab_v3_06_16'));
addpath(genpath('D:/流浪C盘/m_map'));
file='D:/流浪C盘/CMEMS_energy_2014_2023/GLORYS_eddy_fullpath_20160928_20170728_110_123_16_24.nc';
dates=datetime({'2016-10-13','2016-11-17','2016-12-22','2017-01-26','2017-03-09','2017-04-21','2017-06-02','2017-07-14'},'InputFormat','yyyy-MM-dd');
for name={'uo','vo','so','thetao'}
 inf=ncinfo(file,name{1}); assert(isequal({inf.Dimensions.Name},{'longitude','latitude','depth','time'}));
end
results=cell(8,1);
for k=1:8
 fprintf('\nSTAGE %d/8: %s (30-day window)\n',k,string(dates(k)));
 results{k}=stage_kernel(file,dates(k));
 stage=results{k}; save(fullfile(root,sprintf('stage_%s.mat',datestr(dates(k),'yyyy-mm-dd'))),'-struct','stage','-v7.3');
end

%% Export depth-specific panels using a common range across stages; the manuscript uses the 92 m group.
labels={'0m','92m','318m','upper150m'};
for component={'BT','BC'}
 key=component{1};
 for layer=1:4
  all=[];
  for k=1:8, a=results{k}.(key)(:,:,layer); all=[all; a(isfinite(a))]; end
  if isempty(all), continue; end
  lim=prctile(abs(all),98); assert(lim>0);
  folder=fullfile(root,['C' key],labels{layer}); if ~isfolder(folder), mkdir(folder); end
  for k=1:8
   a=results{k}.(key)(:,:,layer);
   if ~any(isfinite(a),'all'), fprintf('SKIPPED invalid layer: %s %s %s\n',key,labels{layer},string(dates(k))); continue; end
   draw_stage(a,results{k}.lon,results{k}.lat,lim,fullfile(folder,sprintf('C%s_%02d_%s_%s.png',key,k,datestr(dates(k),'yyyy-mm-dd'),labels{layer})));
  end
  fprintf('EXPORTED C%s %s: shared limits +/- %.8g W/m^3\n',key,labels{layer},lim);
 end
end
T=table(dates(:),dates(:)-days(15),dates(:)+days(14),'VariableNames',{'TargetDate','WindowStart','WindowEnd'});
writetable(T,fullfile(root,'date_windows.csv'));
fprintf('\nCOMPLETE: MATLAB PNG export finished; invalid N2 layers skipped. Main comparison: CBT/upper150m and CBC/upper150m.\n');
end
function S=stage_kernel(file,target)
% Return conversion fields and coordinates for one target window.
velocityFile=file; hydroFile=file;
lonRange = [110 123];
latRange = [16 23.5];
rho0 = 1025;
g = 9.81;
earthRadius = 6371000;
blockSize = 20;

lonAll = double(ncread(velocityFile,'longitude'));
latAll = double(ncread(velocityFile,'latitude'));
depth = double(ncread(velocityFile,'depth'));
time = ncread(velocityFile,'time');
times=cf_time(velocityFile);
ids=find(times>=target-days(15) & times<=target+days(14));
assert(numel(ids)==30 && all(diff(ids)==1),'Incomplete 30-day window');
firstIndex=ids(1);
iLon = find(lonAll >= lonRange(1) & lonAll <= lonRange(2));
iLat = find(latAll >= latRange(1) & latAll <= latRange(2));
lon = lonAll(iLon); lat = latAll(iLat);
nx = numel(lon); ny = numel(lat); nz = numel(depth); nt = numel(ids);

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
    start = [iLon(1) iLat(1) 1 t0+firstIndex-1];
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

%% Estimate reference stratification within this window.
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
    start = [iLon(1) iLat(1) 1 t0+firstIndex-1];
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


idx=zeros(1,3); for k=1:3, targets=[0 100 300]; [~,idx(k)]=min(abs(depth-targets(k))); end
upper=find(depth<=150); w=layer_weights(depth(upper));
S.BT=cat(3,C_BT(:,:,idx),weighted_vertical_mean(C_BT(:,:,upper),w));
S.BC=cat(3,C_BC(:,:,idx),weighted_vertical_mean(C_BC(:,:,upper),w));
S.lon=lon; S.lat=lat; S.depth=depth; S.N2=N2; S.target=target;
S.actualDepths=depth(idx); S.window=[target-days(15),target+days(14)];
S.C_BT=C_BT; S.C_BC=C_BC;
assert(any(isfinite(S.BT(:,:,4)),'all'));
assert(any(isfinite(S.BC(:,:,4)),'all'));
end
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


function draw_stage(a,lon,lat,lim,filename)
% Export one conversion map with the eddy position.
fig=figure('Visible','off','Color','w','Position',[100 80 820 680]);
m_proj('miller','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
m_pcolor(lon,lat,a'); shading interp; hold on;
caxis([-lim lim]); colormap(jet(256));
cb=colorbar; ylabel(cb,'W m^{-3}');
m_coast('patch',[.8 .8 .8],'edgecolor',[.3 .3 .3]);
m_grid('fancy','tickdir','out');
set(gca,'FontName','Times New Roman','FontSize',10);
tok=regexp(filename,'(20[0-9]{2}-[0-9]{2}-[0-9]{2})','tokens','once'); mark_eddy(datetime(tok{1},'InputFormat','yyyy-MM-dd'));
exportgraphics(fig,filename,'Resolution',300); close(fig);
end

function t=cf_time(file)
% Decode the time coordinate using its NetCDF epoch and units.
unit=char(ncreadatt(file,'time','units')); raw=double(ncread(file,'time'));
tok=regexp(unit,'(milliseconds|seconds|hours|days) since ([0-9]{4}-[0-9]{2}-[0-9]{2})(?:[ T]([0-9]{2}:[0-9]{2}:[0-9]{2}))?','tokens','once');
assert(~isempty(tok),'Unsupported time units'); clock='00:00:00'; if numel(tok)>2 && ~isempty(tok{3}), clock=tok{3}; end
base=datetime([tok{2} ' ' clock],'InputFormat','yyyy-MM-dd HH:mm:ss');
switch tok{1}, case 'milliseconds', t=base+milliseconds(raw); case 'seconds', t=base+seconds(raw); case 'hours', t=base+hours(raw); case 'days', t=base+days(raw); end
t=t(:);
end


function mark_eddy(target)
% Interpolate the tracked eddy position to the target date.
T=readtable(fullfile(fileparts(mfilename('fullpath')),'typical_eddy_track.csv'));
d=days(target-datetime(2016,9,1))+1;
x=interp1(T.day,T.longitude,d,'linear');y=interp1(T.day,T.latitude,d,'linear');
if isfinite(x)&&isfinite(y)
 m_plot(x,y,'ko','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5);
end
end
