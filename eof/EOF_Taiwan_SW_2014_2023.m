% ====================================================================
% Weighted EOF decomposition of daily sea level anomalies southwest of Taiwan.
% Select ocean cells, remove the time mean at each cell, weight by area, and decompose.
% Loadings are in cm per PC standard deviation. The 90-day mean is for display only.
% ====================================================================

% Local EOF analysis of daily SLA southwest of Taiwan, 2014--2023.
% Each exported PNG is a single plot without MATLAB title/subplot.
clear; clc; close all;

infile = 'D:\流浪C盘\sla_total.nc';
outdir = fullfile(getenv('USERPROFILE'),'Desktop', ...
    '图片-去除subplot','EOF_Taiwan_SW');
if ~exist(outdir,'dir'); mkdir(outdir); end

% The eddy catalogue selects genesis near 115--122 E, 20--23.5 N.
% Extend the EOF box one degree south to include adjacent circulation.
box = [115 121.938; 19 23.5]; % [west east; south north], degrees
nmode = 6;
lon0 = double(ncread(infile,'longitude'));
lat0 = double(ncread(infile,'latitude'));
traw = double(ncread(infile,'time'));
date = datetime(1970,1,1) + seconds(traw(:));
assert(numel(date)==3652 && date(1)==datetime(2014,1,1) && ...
    date(end)==datetime(2023,12,31),'Unexpected SLA time coverage');

[lon,lat,Y,water] = read_box(infile,lon0,lat0,box,numel(date));
[map,pc,expvar] = eof_box(Y,lat,water,nmode);
fprintf('Selected box: %.3f--%.3f E, %.3f--%.3f N; %d/%d wet cells\n', ...
    min(lon),max(lon),min(lat),max(lat),sum(water(:)),numel(water));
fprintf('Explained variance (%%), modes 1--6: %s\n', ...
    sprintf('%.4f ',expvar));
fprintf('Cumulative first 2: %.4f%%; first 6: %.4f%%\n', ...
    sum(expvar(1:2)),sum(expvar));

% Check whether slight changes in the western and southern boundaries
% preserve the leading temporal patterns. Sign is arbitrary, hence |r|.
altBoxes = cat(3,[116 121.938;19 23.5], ...
                 [115 121.938;19.5 23.5]);
sens = zeros(2,2,2);
for ib = 1:2
    [~,la2,Y2,wet2] = read_box(infile,lon0,lat0,altBoxes(:,:,ib),numel(date));
    [~,pc2,ev2] = eof_box(Y2,la2,wet2,2);
    R = corr(pc(:,1:2),pc2(:,1:2));
    sens(:,:,ib) = abs(R);
    fprintf('Nearby box %d [%g %g; %g %g]: var1 %.2f, var2 %.2f, |PC correlations| [%s]\n', ...
        ib,altBoxes(1,1,ib),altBoxes(1,2,ib), ...
        altBoxes(2,1,ib),altBoxes(2,2,ib),ev2(1),ev2(2), ...
        sprintf('%.2f ',abs(R(:))));
end

% Save concise reproducibility data and images.
T = table((1:nmode)',expvar(:),'VariableNames',{'Mode','VariancePercent'});
writetable(T,fullfile(outdir,'EOF_variance.csv'));
save(fullfile(outdir,'EOF_Taiwan_SW_results.mat'), ...
    'lon','lat','date','box','map','pc','expvar','sens');

% Blue-white-red diverging colours, appropriate for signed SLA loading.
n = 128;
blue = [linspace(.08,1,n)' linspace(.30,1,n)' ones(n,1)];
red = [ones(n,1) linspace(1,.17,n)' linspace(1,.12,n)'];
cmap = [blue;red];
haveCoast = exist('coastlines.mat','file')==2;
if haveCoast
    c = load('coastlines.mat');
end
for k = 1:nmode
    Z = map(:,:,k)';
    lim = max(abs(Z(isfinite(Z))));
    f = figure('Visible','off','Color','w','Position',[70 70 860 620]);
    ax = axes(f);
    im = imagesc(ax,lon,lat,Z);
    im.AlphaData = double(isfinite(Z));
    set(ax,'YDir','normal','Color',[.84 .84 .84], ...
        'FontName','Times New Roman','FontSize',12,'LineWidth',1);
    axis(ax,'tight'); axis(ax,'equal');
    xlim(ax,[box(1,1) box(1,2)]); ylim(ax,[box(2,1) box(2,2)]);
    colormap(ax,cmap); clim(ax,[-lim lim]);
    cb = colorbar(ax); cb.Label.String = 'SLA (cm per 1 s.d. PC)';
    hold(ax,'on');
    if haveCoast
        plot(ax,c.coastlon,c.coastlat,'k','LineWidth',.8);
    end
    xlabel(ax,'Longitude (°E)'); ylabel(ax,'Latitude (°N)');
    grid(ax,'on'); set(ax,'Box','on');
    exportgraphics(f,fullfile(outdir,sprintf('EOF_SW_mode%02d_map.png',k)), ...
        'Resolution',300);
    close(f);

    f = figure('Visible','off','Color','w','Position',[70 70 1100 440]);
    ax = axes(f);
    plot(ax,date,pc(:,k),'Color',[.62 .69 .78],'LineWidth',.65);
    hold(ax,'on');
    plot(ax,date,movmean(pc(:,k),90,'Endpoints','shrink'), ...
        'Color',[.12 .31 .65],'LineWidth',1.65);
    yline(ax,0,'k-','LineWidth',.7);
    xlim(ax,[date(1) date(end)]);
    xticks(ax,datetime(2014:2:2024,1,1));
    xticklabels(ax,cellstr(string(2014:2:2024)));
    xlabel(ax,'Year'); ylabel(ax,sprintf('PC%d (standard deviations)',k));
    legend(ax,{'Daily','90-day moving mean'},'Location','northeast', ...
        'Box','off');
    set(ax,'FontName','Times New Roman','FontSize',11,'LineWidth',1);
    grid(ax,'on'); set(ax,'Box','on');
    exportgraphics(f,fullfile(outdir,sprintf('EOF_SW_mode%02d_pc.png',k)), ...
        'Resolution',300);
    close(f);
end
fprintf('Saved 12 PNGs and result files in %s\n',outdir);

function [lon,lat,Y,wet] = read_box(file,lon0,lat0,box,nt)
% Extract the regional SLA matrix and its valid ocean-cell mask.
    ix = find(lon0>=box(1,1) & lon0<=box(1,2));
    iy = find(lat0>=box(2,1) & lat0<=box(2,2));
    lon = lon0(ix); lat = lat0(iy);
    A = double(ncread(file,'sla',[ix(1) iy(1) 1], ...
        [numel(ix) numel(iy) nt]));
    A = 100*A; % m to cm
    wet = all(isfinite(A),3);
    Y = reshape(A,[],nt);
    Y = Y(wet(:),:);
    Y = Y - mean(Y,2);
end

function [map,pc,expvar] = eof_box(Y,lat,wet,nmode)
% Y has one row per ocean cell and one column per day. U therefore contains spatial vectors and V temporal vectors.
    [nx,ny] = size(wet);
    w = repmat(sqrt(cosd(lat(:)')),nx,1);
    w = w(wet(:));
    X = Y.*w;
    total = sum(X(:).^2);
    [U,S,V] = svds(X,nmode);
    [ss,order] = sort(diag(S),'descend');
    U = U(:,order); V = V(:,order);
    expvar = 100*(ss.^2)/total;
    % Multiply temporal singular vectors by singular values to obtain PC amplitudes.
    pc0 = V.*ss';
    psd = std(pc0,0,1);
    pc = pc0./psd;
    % Undo sqrt(cos(latitude)) weighting and absorb the PC standard deviation into the loading.
    loading = (U./w).*psd;
    map = nan(nx,ny,nmode);
    for k = 1:nmode
        z = nan(nx*ny,1); z(wet(:)) = loading(:,k);
        z = reshape(z,nx,ny);
        [~,imax] = max(abs(z(:)));
        if z(imax)<0
            z = -z; pc(:,k) = -pc(:,k);
        end
        map(:,:,k) = z;
    end
end
