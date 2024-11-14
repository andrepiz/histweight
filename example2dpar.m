%% Example

addpath('permn')

rng(10)
c = colormap('parula');
% scenario = 'square uniform';
% scenario = 'circle random';
scenario = 'points random';

method = 'area'; % 'area','diff','invsquared','gaussian'
gra = 1; % granularity. Default is 1
threads = 6; % number of threads for parallel computation

R = 15;
v1 = 1;
v2 = 10;
xshift = 20;
yshift = -40;

switch scenario
    case 'square uniform'
        n = 5e4;
        % Uniform points on square with larger intensity at the border
        nside = round(sqrt(n));
        [xcoord, ycoord] = meshgrid(R*linspace(-1, 1, nside), R*linspace(-1, 1, nside));
        xcoord = xcoord(:)' + xshift;
        ycoord = ycoord(:)' + yshift;
        coords = [xcoord; ycoord];
        values = v1*ones(1, length(xcoord));
        % Increase intensity of points at the sides
        values(xcoord - xshift > 0.5*R | ycoord - yshift > 0.5*R) = v2;
    case 'circle random'
        n = 5e4;
        % Random points on circle with larger intensity at the border
        radius = R*sqrt(rand(1, n));
        theta = 2*pi*rand(1, n);
        xcoord = xshift + radius.*cos(theta);
        ycoord = yshift + radius.*sin(theta);
        coords = [xcoord; ycoord];
        values = v1*ones(1, n);
        % Increase intensity of points at the center
        values((xcoord-xshift).^2+(ycoord-yshift).^2<(R/2)^2) = v2;
    case 'points random'
        n = 1e6;
        % 2D scattered points
        xcoord = xshift + R*rand(1, n);
        ycoord = yshift + R*rand(1, n);
        coords = [xcoord; ycoord];
        values = [v1*ones(1, floor(n/2)), v2*ones(1, n-floor(n/2))];
end

limits = [floor(min(coords, [], 2)), 1 + ceil(max(coords, [], 2))];

%%---
tic
[bins_hw, counts_hw, edges_hw] = histweight(coords, values, limits, gra, 'method', method);
t_elapsed_hw = toc
%%--

%%---PARALLEL COMPUTATION
p = gcp('nocreate');
if isempty(p)
    parpool
end
tic
[bins_parhw, counts_parhw, edges_parhw] = parhistweight(coords, values, limits, gra, 'method', method, 'nthreads',threads);
t_elapsed_parhw = toc
%%--

%% PLOT

xbincoords = edges_hw{1};
xbincoords = xbincoords(1:end-1) + 0.5;
ybincoords = edges_hw{2};
ybincoords = ybincoords(1:end-1) + 0.5;
[XX, YY] = ndgrid(xbincoords, ybincoords);

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]), 
colormap(c)

ax1 = subplot(1,3,1);
grid on, hold on, axis equal
xlim(limits(1,:) + [-R/2,R/2])
ylim(limits(2,:) + [-R/2,R/2])
scatter(coords(1,:), coords(2,:), [], values,'o','filled')
c = colorbar;
c.Label.String = 'intensity';
title('sampled points')
xlabel('x')
ylabel('y')

ax2 = subplot(1,3,2);
grid on, hold on, axis equal
xlim(gra*xlim(ax1))
ylim(gra*ylim(ax1))
h1 = imagesc([xbincoords(1), xbincoords(end)], [ybincoords(1), ybincoords(end)], bins_hw');
set(h1, 'AlphaData', bins_hw'~=0)
c = colorbar;
c.Label.String = 'intensity';
title('histweight')
xlabel('x')
ylabel('y')

ax3 = subplot(1,3,3);
grid on, hold on, axis equal
xlim(gra*xlim(ax1))
ylim(gra*ylim(ax1))
h2 = imagesc([xbincoords(1), xbincoords(end)], [ybincoords(1), ybincoords(end)], bins_parhw');
set(h2, 'AlphaData', bins_parhw'~=0)
c = colorbar;
c.Label.String = 'counts';
title('histweight parallelized')
xlabel('x')
ylabel('y')
