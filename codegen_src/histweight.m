function [bins, counts, edges] = histweight(dCoords, dValues, dLimits, dGranularity, i32Method, bFlagProgress, bVECTORIZED, bDEBUG_MODE)%#codegen
% HISTWEIGHT weights and bin scattered data points into uniform quantiles of 
% specified granularity within the specified limits. 
% Each data point is expressed in D-dimensional coordinates and has an 
% associated 1-dimensional value. The dimensions are defined with respect to a 
% D-dimensional grid (D1 = rows, D2 = columns, D3 = pages).
% The limits can be different for each dimension. 
% The granularity downsample the limits and increase the number of quantiles. 
% Each value is spread to the neighbouring 
% bins with a weight defined by three different methods:
%   invsquared: inverse squared distance with each vertex
%   diff: 1 minus distance normalized over maximum distance
%   area: fraction of square box area going to each sector
%
% INPUTS:
%   coords [D x N]
%   values [1 x N]
%   limits [D x 2]
%   granularity [1]
%   method [char]
%
% OUTPUTS:
%   bins [M1 x ... x Mi x ... MD]
%   counts [M1 x ... x Mi x ... MD]
%   edges [M1 x ... x Mi x ... MD]

arguments
    dCoords        (2, :) double {ismatrix}
    dValues        (1, :) double {isvector}
    dLimits        (:, 2) double {isvector} = [floor(min(dCoords,[],2)), 1 + ceil(max(dCoords,[],2))];
    dGranularity   (1, 1) double {isscalar} = 1 % uint32?
    i32Method     (1, :) string            = 'area'
    bFlagProgress  (1, 1) logical           = false;
    bVECTORIZED    (1, 1) logical           = false;
    bDEBUG_MODE    (1, 1) logical           = false;
end

bFlagProgress = true;
bVECTORIZED = false;

% Preliminary checks
ui32CoordRowSize = uint32(size(dCoords, 1));
ui32CoordColSize = uint32(size(dCoords, 2));

assert(dGranularity >= 1, 'Please provide granularity as a scalar integer larger or equal to 1')
assert(size(dValues,1) == 1 && size(dValues,2) == ui32CoordColSize, 'Please provide values as [1xN] vector, where N is the second dimension of coords')
assert(size(dLimits,1) == ui32CoordRowSize && size(dLimits,2) == 2, 'Please provide limits as [Dx2] vector, where D is the first dimension of coords')


% DEVNOTE: coder.ignoreConst
% assert(strcmpi(charMethod, 'area') || strcmpi(charMethod, 'diff') || ...
%     strcmpi(charMethod, 'invsquared'), 'Not Supported. Use area, diff or invsquared methods')

% Shift in case of negative data
dShiftVector = min(0, dLimits(:, 1) - 1);

% Find coordinates of points and centers
dCoordRanges = (dLimits - dShiftVector) * dGranularity; % Find allowed ranges of values
dGridPoints  = (dCoords - dShiftVector) * dGranularity;  % Compute grid points

dCenterPoints = ( 0.5 * sign( dGridPoints - round(dGridPoints) ) + round(dGridPoints)); % Compute centres
dVectorCenters2points = dGridPoints - dCenterPoints;

% Relative difference with neighbouring bins
centers2vertexes = permn([-1, 0, 1], ui32CoordRowSize);
% DEVNOTE: Isn't the output of this function constant if number of rows of coords is fixed?
% DEVNOTE 2: is ui32CoordRowSize fixed for each call of this function? If so, move to caller and pass as
% input or mark as persistent.

% Init
M = dCoordRanges(:, 2); % M seems a duplicated variable: really needed?

if isscalar(M)
    bins   = zeros(M, 1);
    counts = zeros(M, 1, 'uint32');
else
    bins   = zeros(M');
    counts = zeros(M', 'uint32');
end

% bins_check = bins;
% counts_check = counts;


% Define indices and loop sizes
ui32S = uint32(size(centers2vertexes, 1)); % what the hell is this S? I do not get
% dw = zeros(1, ui32S); % what is this?
i32Method = 'area';

dCONST_METHOD_DIFF = 1.5.*sqrt(2); % Constant scaling required by diff method

for ii = 1:ui32CoordColSize
    dw = zeros(1, ui32S); % what is this?

    switch i32Method
        case 'invsquared'
            % Inverse squared distance with each vertex
            % d = vecnorm(centers2vertexes(jj, :)' - dVectorCenters2points(:, ii));
            % dw(jj) = 1./(d.^2);

            % Vectorized code
            d = vecnorm(centers2vertexes' - dVectorCenters2points(:, ii), 2, 1);
            dw = 1./(d.^2);
            dw = dw./sum(dw);

        case 'diff'
            % 1 minus distance normalized over maximum distance
            % d = vecnorm(centers2vertexes(jj, :)' - dVectorCenters2points(:, ii));
            % dw(jj) = 1 - d./dCONST_METHOD_DIFF;

            % Vectorized code
            d = vecnorm(centers2vertexes' - dVectorCenters2points(:, ii), 2, 1);
            dw = 1 - d./dCONST_METHOD_DIFF;
            dw = dw./sum(dw);

        case 'area'

            % Fraction of [1x1] box area going to each sector
            % dSides = centers2vertexes(jj, :)' - dVectorCenters2points(:, ii); % DEVNOTE: there may be some way to vectorize this operation I think

            % if all(abs(dSides) < 1)
            % note: all the weights should sum to 1
            % w(jj) = (|x| + (0.5 - |x|) + (0.5 - |x|))*(|y| + (0.5 - |y|) + (0.5 - |y|))*...
            % dw(jj) = prod(1 - abs(dSides));
            % end % Else what happens?

            % Vectorized mode
            dSides = centers2vertexes' - dVectorCenters2points(:, ii);
            bNonZeroMask = all(abs(dSides) < 1, 1);
            dw(bNonZeroMask) = prod(1 - abs(dSides(:, bNonZeroMask)));

            assert(sum(dw) - 1.00 <= 1.5*eps, 'Area method: weights do not sum to 1')

    end


    % Assigning weighted value to each array element and its neighbours in the D-dimensional space

    % bNonZeroMask % can be reused here to replace loop
    % DEVNOTE: sector_idx may be converted to integers operations. Not sure max min can work with different
    % types though. UPDATE: it can with scalars and perfoms downcasting to integers.
    tic

    if not(bVECTORIZED)
        for jj = 1:ui32S
            if dw(jj) > 0

                % Find index of sector
                sector_idx = ceil(dCenterPoints(:, ii) + centers2vertexes(jj, :)'); % vectorizable: OK

                % Snapping index out of the range to closest neighbours
                sector_idx_sat = max( dCoordRanges(:, 1), min( sector_idx, dCoordRanges(:, 2)) );

                % neighbour_idx_cell = num2cell(sector_idx_sat)';
                idR = sector_idx_sat(1);
                idC = sector_idx_sat(2);

                bins(idR, idC) = bins(idR, idC) + dValues(ii)*dw(jj);

                counts(idR, idC) = counts(idR, idC) + 1;

            end
        end

    else
        % Vectorized version
        centers2vertexes_reduced = centers2vertexes(bNonZeroMask', :);
        dw_reduced = dw(bNonZeroMask);

        % sector_idx_vectorized = zeros(size(centers2vertexes_reduced'), 'double');
        sector_idx_vectorized = ceil(dCenterPoints(:, ii) + centers2vertexes_reduced'); % vectorizable

        % Enforce bounds for sector_idx
        % sector_idx_bounder_vectorized = zeros(size(sector_idx_vectorized), 'double');
        sector_idx_bounder_vectorized = max(dCoordRanges(:, 1), min( sector_idx_vectorized, dCoordRanges(:, 2)) );

        % neighbour_idx_cell_vectorized = num2cell(sector_idx_bounder_vectorized)';
        % dWeightedValues = dValues(ii)*dw;

        % Get allocation linear indices
        linearIdx = sub2ind(dCoordRanges(:, 2)', sector_idx_bounder_vectorized(1, :), sector_idx_bounder_vectorized(2, :));
        
        % Allocate values
        bins(linearIdx) = bins(linearIdx) + dValues(ii).*dw_reduced; % Allocate dWeightedValues for each ii in all indexed entries of bins
        counts(linearIdx) = counts(linearIdx) + uint32(1);

        if bDEBUG_MODE
            % assert(any( (bins - bins_check) <= 1.5*eps, 'all') )
            % assert(all( (counts - counts_check) == 0, 'all') )
        end

    end
    % Print computation progress
    if mod(ii, round(ui32CoordColSize/100)) == 0 && bFlagProgress
        fprintf('     %d%%\n', int32(round(1e2*ii/ui32CoordColSize)));
    end
end

wall_time = toc;
fprintf('\n\tWall time of histweight call %5.4f [ms]', 1000*wall_time);

%   bins [M1 x ... x Mi x ... MD]
%   counts [M1 x ... x Mi x ... MD]
%   edges [M1 x ... x Mi x ... MD]
 
edges = cell(1, ui32CoordRowSize); % Not clear why this should be a cell? DEVNOTE: is that it can have different sized arrays within it?
for kk = 1:ui32CoordRowSize
    edges{kk} = (0:1:size(bins, kk)) + dShiftVector(kk)*dGranularity;
end


end

