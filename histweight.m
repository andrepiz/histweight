function [bins, counts, edges] = histweight(coords, values, limits, granularity, method)
% HISTWEIGHT weights and bin scattered data points into uniform quantiles of 
% specified granularity within the specified limits. 
% Each data point is expressed in D-dimensional coordinates and has an 
% associated 1-dimensional value. The limits can be different for each dimension. 
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
%   granularity [1]
%   limits [D x 2]
%   weightmethod [char]
%
% OUTPUTS:
%   bins [M1 x ... x Mi x ... MD]
%   counts [M1 x ... x Mi x ... MD]
%   edges [M1 x ... x Mi x ... MD]

% Preliminary checks
[D, N] = size(coords);
if nargin == 2
    limits = [floor(min(coords,[],2)), 1 + ceil(max(coords,[],2))];
    granularity = 1;
    method = 'area';
elseif nargin == 3
    granularity = 1;
    method = 'area';
elseif nargin == 4
    method = 'area';
elseif nargin ~= 5
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
    error('Please provide values as [Dx2] vector, where D is the first dimension of coords')
end

% Shift in case of negative data
shift = min(0, limits(:, 1) - 1);

% Find coordinates of points and centers
ranges = (limits - shift)*granularity;
points = (coords - shift)*granularity;
centers = (0.5*sign(points - round(points)) + round(points));
centers2points = points - centers;

% Relative difference with neighbouring bins
centers2vertexes = permn([-1, 0, 1], D);

% Init
M = ranges(:, 2);
if length(M) == 1
    bins = zeros(M, 1);
    counts = zeros(M, 1);
else
    bins = zeros(M');
    counts = zeros(M');
end
disp('   Binning...')
for ii = 1:N
    
    S = size(centers2vertexes, 1);
    w = zeros(S, 1);

    for jj = 1:S
        switch method
            case 'invsquared'
                % Inverse squared distance with each vertex
                d = vecnorm(centers2vertexes(jj, :)' - centers2points(:, ii));
                w(jj) = 1/(d.^2);
            case 'diff'
                % 1 minus distance normalized over maximum distance
                d = vecnorm(centers2vertexes(jj, :)' - centers2points(:, ii));
                w(jj) = 1 - d/(1.5*sqrt(2));
            case 'area'
                % Fraction of [1x1] box area going to each sector
                sides = centers2vertexes(jj, :)' - centers2points(:, ii);
                if all(abs(sides) < 1)
                    % note: all the weights should sum to 1
                    % w(jj) = (|x| + (0.5 - |x|) + (0.5 - |x|))*(|y| + (0.5 - |y|) + (0.5 - |y|))*... 
                    w(jj) = prod(1 - abs(sides));
                end
            otherwise
                error('Not Supported. Use area, diff or invsquared methods')
        end
    end

    if strcmp(method, 'area') 
        if sum(w)-1 > eps
            error('Area method: weights do not sum to 1')
        end
    else
        % normalize weights such that sum to 1
        w = w./sum(w);
    end

    % Assigning weighted value to each array element and its neighbours in the D-dimensional space 
    for jj = 1:S
        if w(jj) > 0 
            % Find index of sector
            sector_idx = ceil(centers(:, ii) + centers2vertexes(jj, :)');
            % Snapping index out of the range to closest neighbours
            sector_idx_sat = max(ranges(:, 1), min(sector_idx, ranges(:, 2)));
            neighbour_idx_cell = num2cell(sector_idx_sat)';
            bins(neighbour_idx_cell{:}) = bins(neighbour_idx_cell{:}) + values(ii)*w(jj);
            counts(neighbour_idx_cell{:}) = counts(neighbour_idx_cell{:}) + 1;
        end
    end
    if mod(ii, round(N/10)) == 0
        disp(['     ',num2str(round(1e2*ii/N)),'%'])
    end
end

edges = cell(1, D);
for kk = 1:D
    edges{kk} = [0:1:size(bins, kk)] + shift(kk)*granularity;
end

end

