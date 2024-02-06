%% Example

addpath('permn')

rng(10)
n = 1e3;
scenario = 'sine';

R = 2;
v1 = 1;
v2 = 2;
shift = 1;

switch scenario
    case 'sine'
        % Uniform points on sine with larger intensity before half
        coords = shift + linspace(0, 10*pi, n);
        values = v1*R*sin(coords);
        % Increase intensity of points before half
        values(coords - shift < 5*pi) = v2*R*sin(coords(coords - shift < 5*pi));
end

limits = [floor(min(coords, [], 2)), 1 + ceil(max(coords, [], 2))];
gra = 3;
[bins_hw, counts_hw, edges_hw] = histweight(coords, values, limits, gra);

xbincoords = edges_hw{1};
xbincoords = xbincoords(1:end-1) + 0.5;
[XX] = ndgrid(xbincoords);

% histcount comparison
bins_hc = histcounts(coords(1,:)*gra, edges_hw{1});

figure(), 

ax1 = subplot(1,3,1);
grid on, hold on
xlim(limits(1,:) + [-R/2,R/2])
plot(coords, values,'o-')
title('sampled points')
xlabel('x')
ylabel('y')

ax2 = subplot(1,3,2);
grid on, hold on
xlim(gra*xlim(ax1))
bar(xbincoords, bins_hw')
title('histweight')
xlabel('x')
ylabel('y')

ax3 = subplot(1,3,3);
grid on, hold on
xlim(gra*xlim(ax1))
bar(xbincoords, bins_hc')
title('histcounts')
xlabel('x')
ylabel('y')