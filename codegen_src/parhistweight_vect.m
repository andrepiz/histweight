function [bins, counts, edges] = parhistweight(dCoords, dValues, dLimits, dGranularity, charMethod, ui8Numthreads)%#codegen
% PARHISTWEIGHT Parallel version of histweights. Computation is performed
% by splitting coords and values into a specified number of pools
arguments
    dCoords        double {ismatrix}
    dValues        (1, :) {isvector}
    dLimits        (:, 2) {isvector}      = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    dGranularity   (1, 1) {isscalar}      = 1 % uint32?
    charMethod     (1, :) string          = 'area'
    ui8Numthreads  (1,1) uint8 {isscalar} = 4
end

% Preliminary checks
[D, N] = size(dCoords); % TODO: what is D and what is N?, both can are integers.

assert(dGranularity >= 1, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) ~= 1 && size(dValues,2) ~= N, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == D && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')

ixs_pools = ceil(linspace(1, N, ui8Numthreads + 1));
npools = length(ixs_pools) - 1;

% Creating pools
coords_pools = cell(1, npools);
values_pools = cell(1, npools);

for ix = 1:npools
    ixs_temp = ixs_pools(ix):ixs_pools(ix + 1);
    coords_pools{ix} = dCoords(:, ixs_temp);
    values_pools{ix} = dValues(:, ixs_temp);
end

% DEVNOTE: check histweight.m 
parfor ix = 1:npools
    [bins_pools{ix}, counts_pools{ix}] = histweight_codegen(coords_pools{ix}, values_pools{ix}, dLimits, dGranularity, charMethod, false);
end

% DEVNOTE: why is seems to me that there are duplicated operations done in histweight whose output is not used???

bins = sum(cat(D + 1, bins_pools{:}), D + 1);
counts = sum(cat(D + 1, counts_pools{:}), D + 1);
edges = cell(1, D);
shift = min(0, dLimits(:, 1) - 1);

% DEVNOTE: how to replace cells: 1) matrix with dim MAX, MAX and ptrs arrays for each entry (what is now kk).
% Best efficiency: single buffer with indices to extract each kk. (use the command to convert index to subscripts?)
for kk = 1:D
    edges{kk} = [0:1:size(bins, kk)] + shift(kk)*dGranularity;
end

end

