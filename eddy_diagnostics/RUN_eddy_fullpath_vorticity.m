function RUN_eddy_fullpath_vorticity
% ====================================================================
% Upper-ocean vorticity diagnostics at eight eddy-track dates.
% Integrate velocity and the density moment over available levels; use full bathymetric depth H in both gradients.
% Input depth increases downward. The minus sign in chi recovers the upward-positive convention. Both terms have units s^-2.
% ====================================================================

root='C:/Users/86182/Desktop/图片-去除subplot/Eddy_fullpath_20261106';
if ~isfolder(root),mkdir(root);end
diary(fullfile(root,'vorticity_log.txt')); cleanup=onCleanup(@()diary('off'));
addpath(genpath('C:/Users/86182/Desktop/gsw_matlab_v3_06_16'));
addpath(genpath('D:/流浪C盘/m_map'));
file='D:/流浪C盘/CMEMS_energy_2014_2023/GLORYS_eddy_fullpath_20160928_20170728_110_123_16_24.nc';
topo='C:/Users/86182/Desktop/GEBCO_2023_6x6min_Global.nc';
dates=datetime({'2016-10-13','2016-11-17','2016-12-22','2017-01-26','2017-03-09','2017-04-21','2017-06-02','2017-07-14'},'InputFormat','yyyy-MM-dd');
lon=double(ncread(file,'longitude'));
lat=double(ncread(file,'latitude'));
z=double(ncread(file,'depth'));
z=z(z<=266.1);
 nz=numel(z);
 nx=numel(lon);
ny=numel(lat);
t=cf_time(file);
[LON,LAT]=ndgrid(lon,lat);
tlon=double(ncread(topo,'lon'));
tlat=double(ncread(topo,'lat'));
ix=find(tlon>=min(lon)-.2 & tlon<=max(lon)+.2);
iy=find(tlat>=min(lat)-.2 & tlat<=max(lat)+.2);
hh=-double(ncread(topo,'Height',[ix(1) iy(1)],[numel(ix) numel(iy)]));
hh(hh<=0)=NaN;
H=interp2(tlon(ix),tlat(iy),hh',LON,LAT);
 invH=1./H;
[dInvDepthDx,dInvDepthDy]=metric_gradient(invH,lon,lat);
[dPlanetaryPVDx,dPlanetaryPVDy]=metric_gradient(2*7.292115e-5*sind(LAT)./H,lon,lat);
p=zeros(nx,ny,nz);for iz=1:nz,p(:,:,iz)=gsw_p_from_z(-z(iz)*ones(nx,ny),LAT);end
APV=nan(nx,ny,8);
JEBAR=APV;
for k=1:8
 it=find(t==dates(k));
 assert(numel(it)==1,'Target date missing');
 st=[1 1 1 it];
  ct=[nx ny nz 1];
 u=double(ncread(file,'uo',st,ct));
 v=double(ncread(file,'vo',st,ct));
 sp=double(ncread(file,'so',st,ct));
 pt=double(ncread(file,'thetao',st,ct));
 sa=gsw_SA_from_SP(sp,p,repmat(LON,1,1,nz),repmat(LAT,1,1,nz));
 rho=gsw_rho(sa,gsw_CT_from_pt(sa,pt),p);
 transportX=vertical_integral(u,z);
 transportY=vertical_integral(v,z);
 % APV = -M dot grad(f/H); M is the upper-layer transport.
 APV(:,:,k)=-transportX.*dPlanetaryPVDx-transportY.*dPlanetaryPVDy;
 chi=-9.81/1025*vertical_integral((rho-1025).*reshape(z,1,1,nz),z);
 % JEBAR = J(chi,1/H), with J(a,b) = a_x*b_y - a_y*b_x.
 [dChiDx,dChiDy]=metric_gradient(chi,lon,lat);
 JEBAR(:,:,k)=dChiDx.*dInvDepthDy-dChiDy.*dInvDepthDx;
 fprintf('Calculated %s, deepest standard layer %.5f m\n',string(dates(k)),z(end));
end
save(fullfile(root,'vorticity_all.mat'),'APV','JEBAR','H','lon','lat','z','dates','-v7.3');
fields={APV,JEBAR}; names={'APV','JEBAR'};
for q=1:2
 out=fullfile(root,names{q});if ~isfolder(out),mkdir(out);end
 all=fields{q};lim=prctile(abs(all(isfinite(all))),98);assert(lim>0);
 for k=1:8
  f=figure('Visible','off','Color','w','Position',[100 80 820 680]);
  m_proj('miller','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
  m_pcolor(lon,lat,all(:,:,k)');shading interp;hold on;
  colormap(jet(256));caxis([-lim lim]);cb=colorbar;ylabel(cb,'s^{-2}');
  m_coast('patch',[.8 .8 .8],'edgecolor',[.3 .3 .3]);m_grid('fancy','tickdir','out');
  mark_eddy(dates(k));
  exportgraphics(f,fullfile(out,sprintf('%s_%02d_%s.png',names{q},k,datestr(dates(k),'yyyy-mm-dd'))),'Resolution',300);close(f);
 end
end
fprintf('COMPLETE vorticity export\n');
end

function [dFdx,dFdy] = metric_gradient(F,lon,lat)
% Centered horizontal differences use distances in metres; unsupported edges remain NaN.
Re = 6371000;
lon = double(lon(:)); lat = double(lat(:));
dLambda = mean(diff(deg2rad(lon)));
dPhi = mean(diff(deg2rad(lat)));
dx = Re*cosd(lat).'*dLambda;
dy = Re*dPhi;
dFdx = nan(size(F)); dFdy = nan(size(F));
dFdx(2:end-1,:) = (F(3:end,:)-F(1:end-2,:))./(2*dx);
dFdy(:,2:end-1) = (F(:,3:end)-F(:,1:end-2))./(2*dy);
end

function integralValue = vertical_integral(field,depth)
% Trapezoidal integration uses only adjacent pairs with two valid samples.
depth = double(depth(:));
dz = reshape(diff(depth),1,1,[]);
v1 = field(:,:,1:end-1);
v2 = field(:,:,2:end);
validPair = isfinite(v1) & isfinite(v2);
layerValue = 0.5*(v1+v2).*dz;
layerValue(~validPair) = 0;
integralValue = sum(layerValue,3);
integralValue(~any(validPair,3)) = NaN;
end

function t=cf_time(file)
% Decode NetCDF time relative to the epoch in its units attribute.
unit=char(ncreadatt(file,'time','units')); raw=double(ncread(file,'time'));
tok=regexp(unit,'(milliseconds|seconds|hours|days) since ([0-9]{4}-[0-9]{2}-[0-9]{2})(?:[ T]([0-9]{2}:[0-9]{2}:[0-9]{2}))?','tokens','once');
assert(~isempty(tok),'Unsupported time units'); clock='00:00:00'; if numel(tok)>2 && ~isempty(tok{3}), clock=tok{3}; end
base=datetime([tok{2} ' ' clock],'InputFormat','yyyy-MM-dd HH:mm:ss');
switch tok{1}, case 'milliseconds', t=base+milliseconds(raw); case 'seconds', t=base+seconds(raw); case 'hours', t=base+hours(raw); case 'days', t=base+days(raw); end
t=t(:);
end


function mark_eddy(target)
% Interpolate the tracked eddy position to the mapped day.
T=readtable(fullfile(fileparts(mfilename('fullpath')),'typical_eddy_track.csv'));
d=days(target-datetime(2016,9,1))+1;
x=interp1(T.day,T.longitude,d,'linear');y=interp1(T.day,T.latitude,d,'linear');
if isfinite(x)&&isfinite(y)
 m_plot(x,y,'ko','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5);
end
end
