%% Example

addpath('permn')

rng(10)
c = colormap('parula');
%scenario = 'emisphere random';
scenario = 'emisphere uniform';

gra = 1; % granularity. Default is 1

R = 15;
v1 = 1;
xshift = 20;
yshift = -40;
zshift = -10;

switch scenario
    case 'emisphere random'
        n = 2e5;
        % Random points on emisphere with larger intensity as radius is
        % increasing
        xcoord = xshift + R*(1-2*rand(1, n));
        ycoord = yshift + R*(1-2*rand(1, n));
        zcoord = zshift + R*(1-2*rand(1, n));
        ixs_outside = (xcoord-xshift).^2+(ycoord-yshift).^2+(zcoord-zshift).^2 > R^2 | (ycoord-yshift) < 0;
        xcoord(ixs_outside) = [];        
        ycoord(ixs_outside) = [];
        zcoord(ixs_outside) = [];
        coords = [xcoord; ycoord; zcoord];
        values = v1/R*((xcoord-xshift).^2+(ycoord-yshift).^2+(zcoord-zshift).^2);
    case 'emisphere uniform'
        n = 2e5;
        % Uniform points on emisphere with larger intensity as radius is
        % increasing
        xcoord = xshift + R*(linspace(-1,1,floor(n^(1/3))));
        ycoord = yshift + R*(linspace(-1,1,floor(n^(1/3))));
        zcoord = zshift + R*(linspace(-1,1,floor(n^(1/3))));
        [xcoord, ycoord, zcoord] = ndgrid(xcoord, ycoord, zcoord);
        ixs_outside = (xcoord(:)-xshift).^2+(ycoord(:)-yshift).^2+(zcoord(:)-zshift).^2 > R^2 | (ycoord(:)-yshift) < 0;
        xcoord(ixs_outside) = [];        
        ycoord(ixs_outside) = [];
        zcoord(ixs_outside) = [];
        coords = [xcoord; ycoord; zcoord];
        values = v1/R*((xcoord-xshift).^2+(ycoord-yshift).^2+(zcoord-zshift).^2);
end

limits = [floor(min(coords, [], 2)), 1 + ceil(max(coords, [], 2))];

%%---
[bins_hw, counts_hw, edges_hw] = histweight(coords, values, limits, gra);
%%--

% You can also simply call:
%   [bins_hw, counts_hw, edges_hw] = histweight(coords, values);
% Granularity will be set to 1 as default and limits are automatically
% computed

%% PLOT

xbincoords = edges_hw{1};
xbincoords = xbincoords(1:end-1) + 0.5;
ybincoords = edges_hw{2};
ybincoords = ybincoords(1:end-1) + 0.5;
zbincoords = edges_hw{3};
zbincoords = zbincoords(1:end-1) + 0.5;
[XX, YY, ZZ] = ndgrid(xbincoords, ybincoords, zbincoords);

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]), 
colormap(c)

ax1 = subplot(1,2,1);
grid on, hold on, axis equal
xlim(limits(1,:) + [-R/2,R/2])
ylim(limits(2,:) + [-R/2,R/2])
zlim(limits(3,:) + [-R/2,R/2])
scatter3(xcoord, ycoord, zcoord, [], values, 'o')
c = colorbar;
c.Label.String = 'intensity';
title('sampled points')
xlabel('x')
ylabel('y')
zlabel('z')
view([50 10])

ax2 = subplot(1,2,2);
grid on, hold on, axis equal
xlim(gra*xlim(ax1))
ylim(gra*ylim(ax1))
zlim(gra*zlim(ax1))
scatter3(XX(bins_hw~=0), YY(bins_hw~=0), ZZ(bins_hw~=0), 130, bins_hw(bins_hw~=0), 's','filled')
c = colorbar;
c.Label.String = 'intensity';
title('histweight')
xlabel('x')
ylabel('y')
zlabel('z')
view([50 10])

%% ERROR
disp(['Sum of values: ',num2str(sum(values))])
disp(['Sum of bins: ',num2str(sum(bins_hw,'all'))])