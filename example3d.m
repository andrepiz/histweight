%% Example

addpath('permn')

rng(10)
cm_map = 'parula';
%scenario = 'emisphere random';
scenario = 'emisphere uniform';

method = 'area'; % 'area','diff','invsquared','gaussian'
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
        values = v1/R*((xcoord-xshift).^2+(ycoord-yshift).^2+(zcoord-zshift).^2);
end

xyzcoords = [xcoord; ycoord; zcoord];
xyzlimits = [floor(min(xyzcoords, [], 2)), 1 + ceil(max(xyzcoords, [], 2))];
ijkcoords = [ycoord; xcoord; zcoord]; % defined with respect to 3D matrix
ijklimits = [xyzlimits(2,:); xyzlimits(1,:); xyzlimits(3,:)]; % defined with respect to 3D matrix

%%---
[bins_hw, counts_hw, edges_hw] = histweight(ijkcoords, values, ijklimits, gra, method);
%%--

% You can also simply call:
%   [bins_hw, counts_hw, edges_hw] = histweight(ijkcoords, values);
% Granularity will be set to 1 as default and limits are automatically
% computed. Area method is used as default

%% PLOT

ibincoords = edges_hw{1};
ibincoords = ibincoords(1:end-1) + 0.5;
jbincoords = edges_hw{2};
jbincoords = jbincoords(1:end-1) + 0.5;
kbincoords = edges_hw{3};
kbincoords = kbincoords(1:end-1) + 0.5;
[II, JJ, KK] = ndgrid(ibincoords, jbincoords, kbincoords);
% Remove empty voxels
xbincoords = JJ(bins_hw~=0);
ybincoords = II(bins_hw~=0);
zbincoords = KK(bins_hw~=0);
bins_hw_filtered = bins_hw(bins_hw~=0);

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]), 
colormap(cm_map)

ax1 = subplot(1,2,1);
grid on, hold on, axis equal
xlim(xyzlimits(1,:) + [-R/2,R/2])
ylim(xyzlimits(2,:) + [-R/2,R/2])
zlim(xyzlimits(3,:) + [-R/2,R/2])
scatter3(xcoord, ycoord, zcoord, [], values, 'o')
cb = colorbar;
cb.Label.String = 'intensity';
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
scatter3(xbincoords, ybincoords, zbincoords, 130, bins_hw_filtered, 's','filled')
cb = colorbar;
cb.Label.String = 'intensity';
title('histweight')
xlabel('x')
ylabel('y')
zlabel('z')
view([50 10])

%% ERROR
disp(['Sum of values: ',num2str(sum(values))])
disp(['Sum of bins: ',num2str(sum(bins_hw,'all'))])