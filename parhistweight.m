function [bins, counts, edges] = parhistweight(coords, values, limits, granularity, method, nthreads)
% PARHISTWEIGHT Parallel version of histweights. Computation is performed
% by splitting coords and values into a specified number of pools

% Preliminary checks
[D, N] = size(coords);
if nargin == 2
    limits = [floor(min(coords,[],2)), 1 + ceil(max(coords,[],2))];
    granularity = 1;
    method = 'area';
    nthreads = 4;
elseif nargin == 3
    granularity = 1;
    method = 'area';
    nthreads = 4;
elseif nargin == 4
    method = 'area';
    nthreads = 4;
elseif nargin == 5
    nthreads = 4;
elseif nargin ~= 6
    error('Not enough inputs')
end

% Size checks
if length(granularity) ~= 1 || granularity < 1 || floor(granularity) ~= granularity
    error('Please provide granularity as a scalar integer larger or equal to 1')
end
if size(values,1) ~= 1 || size(values,2) ~= N
    error('Please provide values as [1xN] vector, where N is the second dimension of coords')
end
if size(limits,1) ~= D || size(limits,2) ~= 2
    error('Please provide limits as [Dx2] vector, where D is the first dimension of coords')
end

ixs_pools = ceil(linspace(1, N, nthreads + 1));
npools = length(ixs_pools) - 1;
% Creating pools
coords_pools = cell(1, npools);
values_pools = cell(1, npools);
for ix = 1:npools
    ixs_temp = ixs_pools(ix):ixs_pools(ix + 1);
    coords_pools{ix} = coords(:, ixs_temp);
    values_pools{ix} = values(:, ixs_temp);
end

parfor ix = 1:npools
    [bins_pools{ix}, counts_pools{ix}] = histweight(coords_pools{ix}, values_pools{ix}, limits, granularity, method, false);
end

bins = sum(cat(D + 1, bins_pools{:}), D + 1);
counts = sum(cat(D + 1, counts_pools{:}), D + 1);
edges = cell(1, D);
shift = min(0, limits(:, 1) - 1);
for kk = 1:D
    edges{kk} = [0:1:size(bins, kk)] + shift(kk)*granularity;
end
end

