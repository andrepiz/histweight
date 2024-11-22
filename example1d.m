%% Example

addpath('permn')

rng(10)
scenario = 'sine uniform';

method = 'gaussian'; % 'area','diff','invsquared','gaussian'
gra = 1; % granularity. Default is 1

R = 2;
v1 = 1;
v2 = 2;
shift = 1;

switch scenario
    case 'sine uniform'
        n = 1e3;
        % Uniform points on sine with larger intensity before half
        coords = shift + linspace(0, 10*pi, n);
        values = v1*R*sin(coords);
        % Increase intensity of points before half
        values(coords - shift < 5*pi) = v2*R*sin(coords(coords - shift < 5*pi));
end

limits = [floor(min(coords, [], 2)), 1 + ceil(max(coords, [], 2))];

%%---
tic
[bins_hw_GT, counts_hw_GT, edges_hw_GT] = histweight(coords, values, limits, gra, 'method', method);
toc
%%--

methodID = int32(2); % 'area'
bFlagProgress = false;
bVECTORIZED   = false;
bDEBUG_MODE   = false;

addpath("codegen_src/")

tic
[bins_hw_vect, counts_hw_vect, edges_hw_vect] = histweight_vect(coords, values, limits, gra, methodID, bFlagProgress, bVECTORIZED, bDEBUG_MODE);
toc

return

% You can also simply call:
%   [bins_hw, counts_hw, edges_hw] = histweight(coords, values);
% Granularity will be set to 1 as default and limits are automatically
% computed. Area method is used as default.

% histcount comparison
bins_hc = histcounts(coords(1,:)*gra, edges_hw{1});

%% PLOT

xbincoords = edges_hw{1};
xbincoords = xbincoords(1:end-1) + 0.5;
[XX] = ndgrid(xbincoords);

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]), 

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
bar(xbincoords, bins_hw)
title('histweight')
xlabel('x')
ylabel('y')

ax3 = subplot(1,3,3);
grid on, hold on
xlim(gra*xlim(ax1))
bar(xbincoords, bins_hc)
title('histcounts')
xlabel('x')
ylabel('y')

%% ERROR
disp(['Sum of values: ',num2str(sum(values))])
disp(['Sum of bins: ',num2str(sum(bins_hw,'all'))])
