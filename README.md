The MATLAB function histcounts allows to bin scattered data points into quantiles and counts the number of points falling into each quantile.
HISTWEIGHT improves histcounts by considering intensity values for each data point that are spread across the neighbouring quantiles. 
Here's an example of comparison using histweight versus histcounts with randomly scattered data points with the same intensity:
![random_points_different_intensity](https://github.com/andrepiz/histweight/assets/75851004/aaac8b2b-35c5-45cd-9840-3d4819e2f0e8)
And using data points with different intensity:
![random_points_same_intensity](https://github.com/andrepiz/histweight/assets/75851004/88a5992e-4b46-4b29-b31b-4c2d46506751)


**[bins, counts, edges] = histweight(coords, values, limits, granularity, method)**

HISTWEIGHT weights and bin scattered data points into uniform quantiles of 
specified granularity within the specified limits. 
Each data point is expressed in D-dimensional coordinates and has an 
associated 1-dimensional value. The limits can be different for each dimension. 
The granularity downsample the limits and increase the number of quantiles. 
Each value is spread to the neighbouring 
bins with a weight defined by three different methods:
- invsquared: inverse squared distance with each vertex
- diff: 1 minus distance normalized over maximum distance
- area: fraction of square box area going to each sector

**INPUTS**:
- coords [D x N]
- values [1 x N]
- granularity [1]
- limits [D x 2]
- weightmethod [char]

**OUTPUTS**:
- bins [M1 x ... x Mi x ... MD]
- counts [M1 x ... x Mi x ... MD]
- edges [M1 x ... x Mi x ... MD]

