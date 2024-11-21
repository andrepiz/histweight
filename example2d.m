%% Example
clear
clc
close all

%% OPTIONS
bRUN_BENCHMARK = true;
ui32Ntrials = 100;
bPROFILE_CODE = false;

addpath('permn')
addpath("codegen_src/")

%% EXAMPLE SETTINGS
rng(10)
cm_map = 'parula';
%scenario = 'square uniform';
%scenario = 'circle random';
scenario = 'points random';

method = 'area'; % 'area','diff','invsquared','gaussian'
size_kernel = 1;
gra = 1; % granularity. Default is 1
nthreads = 4;
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
        values = v1*ones(1, n);
        % Increase intensity of points at the center
        values((xcoord-xshift).^2+(ycoord-yshift).^2<(R/2)^2) = v2;
    case 'points random'
        n = 1e5;
        % 2D scattered points
        xcoord = xshift + R*rand(1, n);
        ycoord = yshift + R*rand(1, n);
        values = [v1*ones(1, floor(n/2)), v2*ones(1, n-floor(n/2))];
end

% Define 2D inputs
xycoords = [xcoord; ycoord];
xylimits = [floor(min(xycoords, [], 2)), 1 + ceil(max(xycoords, [], 2))];
ijcoords = [ycoord; xcoord]; % defined with respect to 2D matrix
ijlimits = [xylimits(2,:); xylimits(1,:)]; % defined with respect to 2D matrix

% Instantiate arrays for timings
timings_original = nan(1, ui32Ntrials);
timings_original_parallel = nan(1, ui32Ntrials);
timings_optimized = nan(1, ui32Ntrials);
timings_optimized_parallel = nan(1, ui32Ntrials);
timings_optimized_vect = nan(1, ui32Ntrials);
timings_optimized_mex = nan(1, ui32Ntrials);
timings_optimized_mex_vect = nan(1, ui32Ntrials);


% Inputs for histweight_2d and variants
methodID      = int32(2); % 'area'
bFlagProgress = false;
bVECTORIZED   = false;
bDEBUG_MODE   = false;


% DEFINE for MEX versions (required, no optionals due to coder limitation)
dGaussianSigma = 1/3;
dWindowSize = 1;

%% ORIGINAL VERSION
%%---

if bPROFILE_CODE == true
    profStatus = profile('status');
    if strcmpi(profStatus.ProfilerStatus, 'on')
        profile off
    end

    profile on -history -timer performance
end

for idTrial = 1:ui32Ntrials
    tic
    [bins_hw, counts_hw, edges_hw] = histweight(ijcoords, values, ijlimits, gra, 'method', method, 'window', size_kernel);
    timings_original(idTrial) = toc;
end

%% ORIGINAL PARALLEL
for idTrial = 1:ui32Ntrials
    % tic
    % [bins_hw_par, counts_hw_par, edges_hw_par] = parhistweight(ijcoords, values, ijlimits, gra, 'method', method, 'nthreads',nthreads);
    % timings_original_parallel(idTrial) = toc;
end
%%--

if bPROFILE_CODE == true
    profile off
    profile_original = profile('info');
    profsave(profile_original, 'profile_original')
    profile viewer
end

mean_timing_original = mean(timings_original, 'all', "omitnan");
fprintf('\nMean time, original: %4.4g [ms]\n', mean_timing_original/1000)

mean_timing_original_parallel = mean(timings_original_parallel, 'all', "omitnan");
fprintf('\nMean time, original parallelized: %4.4g [ms]\n', mean_timing_original_parallel/1000)


%% OPTIMIZED VERSION NON-VECT
if bPROFILE_CODE == true
    profile on -history -timer performance 
end

for idTrial = 1:ui32Ntrials

    tic
    [bins_hw_2dnoVect, counts_hw_2dnoVect, edges_hw_2dnoVect] = histweight_2d(ijcoords, values, ijlimits, gra, ...
        methodID, bFlagProgress, false, bDEBUG_MODE, dGaussianSigma, dWindowSize);
    timings_optimized(idTrial) = toc;
end

for idTrial = 1:ui32Ntrials
    % tic
    % [bins_hw_2dnoVect_par, counts_hw_2dnoVect_par, edges_hw_2dnoVect_par] = parhistweight_2d(ijcoords, values, ...
    %     ijlimits, gra, methodID, nthreads, dGaussianSigma, dWindowSize);
    % timings_optimized_parallel(idTrial) = toc;
end
%%--

if bPROFILE_CODE == true
    profile off
    profile_optNonVect = profile('info');
    profsave(profile_optNonVect, 'profile_optimized_nonVect')
end

mean_timings_optimized = mean(timings_optimized, 'all', "omitnan");
fprintf('\nMean time, optimized non-vect: %4.4g [ms]\n', mean_timings_optimized/1000)

mean_timing_optimized_parallel = mean(timings_optimized_parallel, 'all', "omitnan");
fprintf('\nMean time, optimized parallelized: %4.4g [ms]\n', mean_timing_optimized_parallel/1000)

%% OPTIMIZED VERSION VECT
if bRUN_BENCHMARK == true

    for idTrial = 1:ui32Ntrials

        tic
        [bins_hw_2dvect, counts_hw_2dvect, edges_hw_2dvect] = histweight_2d(ijcoords, values, ijlimits, gra, methodID, ...
            bFlagProgress, true, bDEBUG_MODE, dGaussianSigma, dWindowSize);
        timings_optimized_vect(idTrial) = toc;
    end

    mean_timings_optimized_vect = mean(timings_optimized_vect, 'all', "omitnan");
    fprintf('\nMean time, optimized vect: %4.4g [ms]\n', mean_timings_optimized_vect/1000)

end

%%  OPTIMIZED MEX VERSION NON-VECT
for idTrial = 1:ui32Ntrials

    tic
    [bins_hw_2dnoVect_MEX, counts_hw_2dnoVect_MEX, edges_hw_2dnoVect_MEX] = histweight_2d_MEX(ijcoords, values, ijlimits, gra, int8(methodID), ...
        bFlagProgress, false, bDEBUG_MODE, dGaussianSigma, dWindowSize);
    timings_optimized_mex(idTrial) = toc;
end

mean_timings_optimized_mex = mean(timings_optimized_mex, 'all', "omitnan");
fprintf('\nMean time, optimized_mex non-vect: %4.4g [ms]\n', mean_timings_optimized_mex/1000)

%%  OPTIMIZED MEX VERSION VECT

if bRUN_BENCHMARK == true

    for idTrial = 1:ui32Ntrials

        tic
        [bins_hw_2dvect_MEX, counts_hw_2dvect_MEX, edges_hw_2dvect_MEX] = histweight_2d_MEX(ijcoords, values, ijlimits, gra, int8(methodID), ...
            bFlagProgress, true, bDEBUG_MODE, dGaussianSigma, dWindowSize);
        timings_optimized_mex_vect(idTrial) = toc;
    end

    mean_timings_optimized_mex_vect = mean(timings_optimized_mex_vect, 'all', "omitnan");
    fprintf('\nMean time, optimized_mex vect: %4.4g [ms]\n', mean_timings_optimized_mex_vect/1000)

end

% You can also simply call:
%   [bins_hw, counts_hw, edges_hw] = histweight(ijcoords, values);
% Granularity will be set to 1 as default and limits are automatically
% computed. Area method is used by default

% histcount comparison (original version)
bins_hc = histcounts2(ijcoords(1,:)*gra, ijcoords(2,:)*gra, edges_hw{1}, edges_hw{2});

% EQUIVALENCE TEST (2d optimized)

binsError_2doptNoVect = bins_hw_2dnoVect - bins_hw;
binsError_2dmexNoVect = bins_hw_2dnoVect_MEX - bins_hw;

fprintf('\n------------------ EQUIVALENCE TEST -------------------\n')
fprintf('\nOptimized version for 2d inputs (no mex):\n')
fprintf('\tSum error: %4.4g\n', sum(binsError_2doptNoVect, 'all'))
fprintf('\tMax error: %4.4g\n', max(binsError_2doptNoVect, [], 'all'))
fprintf('\tMean error: %4.4g\n', mean(binsError_2doptNoVect, 'all'))

fprintf('\nOptimized mexed version for 2d inputs:\n')
fprintf('\tSum error: %4.4g\n', sum(binsError_2dmexNoVect, 'all'))
fprintf('\tMax error: %4.4g\n', max(binsError_2dmexNoVect, [], 'all'))
fprintf('\tMean error: %4.4g\n', mean(binsError_2dmexNoVect, 'all'))


return
%% PLOT

ibincoords = edges_hw{1};
ibincoords = ibincoords(1:end-1) + 0.5;
jbincoords = edges_hw{2};
jbincoords = jbincoords(1:end-1) + 0.5;
xbincoords = jbincoords;
ybincoords = ibincoords;

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]), 
colormap(cm_map)

ax1 = subplot(1,3,1);
grid on, hold on, axis equal
xlim(xylimits(1,:) + [-R/2,R/2])
ylim(xylimits(2,:) + [-R/2,R/2])
scatter(xycoords(1,:), xycoords(2,:), [], values,'o','filled')
cb = colorbar;
cb.Label.String = 'intensity';
title('sampled points')
xlabel('x')
ylabel('y')

ax2 = subplot(1,3,2);
grid on, hold on, axis equal
xlim(gra*xlim(ax1))
ylim(gra*ylim(ax1))
h1 = imagesc([xbincoords(1), xbincoords(end)], [ybincoords(1), ybincoords(end)], bins_hw);
set(h1, 'AlphaData', bins_hw~=0)
cb = colorbar;
cb.Label.String = 'intensity';
title('histweight')
xlabel('x')
ylabel('y')

ax3 = subplot(1,3,3);
grid on, hold on, axis equal
xlim(gra*xlim(ax1))
ylim(gra*ylim(ax1))
h2 = imagesc([xbincoords(1), xbincoords(end)], [ybincoords(1), ybincoords(end)], bins_hc);
set(h2, 'AlphaData', bins_hc~=0)
cb = colorbar;
cb.Label.String = 'counts';
title('histcounts')
xlabel('x')
ylabel('y')

%% ERROR
disp(['Sum of values: ',num2str(sum(values))])
disp(['Sum of bins: ',num2str(sum(bins_hw,'all'))])
disp(['Error: ',num2str(sum(bins_hw,'all') - sum(values))])
