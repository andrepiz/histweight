function [bins, counts, edges] = parhistweight_2d(dCoords, dValues, dLimits, dGranularity, i32Method, ui8Numthreads, dGaussianSigma, dWindowSize)%#codegen
% PARHISTWEIGHT Parallel version of histweights. Computation is performed
% by splitting coords and values into a specified number of pools
arguments
    dCoords        (2, :) double {ismatrix}
    dValues        (1, :) double {isvector}
    dLimits        (:, 2) double {isvector} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    dGranularity   (1, 1) double {isscalar} = 1 % uint32?
    i32Method      (1, :) int32           = 2 % area
    ui8Numthreads  (1,1) uint8 {isscalar} = 4
    dGaussianSigma (1, 1) double {isscalar} = 1/3
    dWindowSize    (1, 1) double {isscalar} = 1;
end


% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')


idx_Pools = ceil(linspace(1, double(ui32CoordColSize), double(ui8Numthreads + uint8(1))));
ui8NumOfPools = length(idx_Pools) - 1;

% Creating pools
coords_pools = cell(1, ui8NumOfPools);
values_pools = cell(1, ui8NumOfPools);

for idx = 1:ui8NumOfPools
    idx_temp = idx_Pools(idx):idx_Pools(idx + 1);
    coords_pools{idx} = dCoords(:, idx_temp);
    values_pools{idx} = dValues(:, idx_temp);
end

% DEVNOTE: check histweight.m 
parfor idx = 1:ui8NumOfPools
    [bins_pools{idx}, counts_pools{idx}] = histweight_2d(coords_pools{idx}, values_pools{idx}, dLimits, dGranularity, i32Method, ...
                                                       false, false, false, dGaussianSigma, dWindowSize);
end

% DEVNOTE: why is seems to me that there are duplicated operations done in histweight whose output is not used???

bins = sum(cat(ui32CoordRowSize + 1, bins_pools{:}), ui32CoordRowSize + 1);
counts = sum(cat(ui32CoordRowSize + 1, counts_pools{:}), ui32CoordRowSize + 1);
edges = cell(1, ui32CoordRowSize);
shift = min(0, dLimits(:, 1) - 1);

% DEVNOTE: how to replace cells: 1) matrix with dim MAX, MAX and ptrs arrays for each entry (what is now kk).
% Best efficiency: single buffer with indices to extract each kk. (use the command to convert index to subscripts?)
for kk = 1:ui32CoordRowSize
    edges{kk} = (0:1:size(bins, kk)) + shift(kk)*dGranularity;
end

end

