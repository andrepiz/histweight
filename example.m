%% Example

rng(10)
n = 5e4;
c = colormap('parula');
scenario = 'circle';
%scenario = 'square';
%scenario = 'random points';

R = 15;
v1 = 1;
v2 = 10;
xshift = 20;
yshift = -40;

switch scenario
    case 'square'
        % Uniform points on square with larger intensity at the border
        nside = round(sqrt(n));
        [xcoord, ycoord] = meshgrid(R*linspace(-1, 1, nside), R*linspace(-1, 1, nside));
        xcoord = xcoord(:)' + xshift;
        ycoord = ycoord(:)' + yshift;
        coords = [xcoord; ycoord];
        values = v1*ones(1, length(xcoord));
        % Increase intensity of points at the sides
        values(xcoord - xshift > 0.5*R | ycoord - yshift > 0.5*R) = v2;
    case 'circle'
        % Random points on circle with larger intensity at the border
        radius = R*sqrt(rand(1, n));
        theta = 2*pi*rand(1, n);
        xcoord = xshift + radius.*cos(theta);
        ycoord = yshift + radius.*sin(theta);
        coords = [xcoord; ycoord];
        values = v1*ones(1, n);
        % Increase intensity of points at the center
        values((xcoord-xshift).^2+(ycoord-yshift).^2<(R/2)^2) = v2;
    case 'random points'
        % 2D scattered points
        xcoord = xshift + R*rand(1, n);
        ycoord = yshift + R*rand(1, n);
        coords = [xcoord; ycoord];
        values = [v1*ones(1, floor(n/2)), v2*ones(1, n-floor(n/2))];
end

% Normal call
%[bins_hw, counts_hw, edges_hw] = histweight(coords, values);

% Changing limits and/or granularity
limits = [floor(min(coords, [], 2)), 1 + ceil(max(coords, [], 2))];
gra = 3;
[bins_hw, counts_hw, edges_hw] = histweight(coords, values, limits, gra);

xbincoords = edges_hw{1};
xbincoords = xbincoords(1:end-1) + 0.5;
ybincoords = edges_hw{2};
ybincoords = ybincoords(1:end-1) + 0.5;
[XX, YY] = ndgrid(xbincoords, ybincoords);

% histcount comparison
bins_hc = histcounts2(coords(1,:)*gra, coords(2,:)*gra, edges_hw{1}, edges_hw{2});

figure(), 
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
h2 = imagesc([xbincoords(1), xbincoords(end)], [ybincoords(1), ybincoords(end)], bins_hc');
set(h2, 'AlphaData', bins_hc'~=0)
c = colorbar;
c.Label.String = 'counts';
title('histcounts')
xlabel('x')
ylabel('y')

%%
sum_values = sum(values, 'all');
sum_bins_hw = sum(bins_hw, 'all');
sum_bins_hc = sum(bins_hc,'all');

err1 = sum_bins_hw - sum_values
err2 = sum_bins_hc - sum_values