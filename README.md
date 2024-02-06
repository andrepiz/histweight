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
