**HISTWEIGHT** bins scattered data points defined in any dimension weighting them into uniform bins of specified granularity. 
Each data point has an associated intensity value which is spread to the neighbouring bins with a weight defined by three different methods (invsquared, diff or area).
Limits of binning can be defined, different for each dimension. The granularity downsample the limits and increase the number of bins. 

Examples of how the function works are depicted in the following figure. Granularity is set to 1 and area method is used for all of them.
- 1-dimensional uniform points sampled from a sine wave with decreasing amplitude:
  ![sine_uniform](https://github.com/andrepiz/histweight/assets/75851004/56afe971-dcec-43c9-86b6-a13621537e0a)
- 2-dimensional random points with two different intensities:
  ![points_random](https://github.com/andrepiz/histweight/assets/75851004/1e2d00f1-c823-49dc-b04d-40e37b5f7da4)
- 2-dimensional random points sampled within a circle, with a larger intensity in a inner circle:
  ![circle_random](https://github.com/andrepiz/histweight/assets/75851004/b01b1dbc-c12b-4688-8932-7eee7ef04a15)
- 3-dimensional uniform points sampled within a emisphere, with intensity increasing with the radius:
  ![emisphere_uniform](https://github.com/andrepiz/histweight/assets/75851004/77ccf275-854f-4b1d-aa4e-88adcade1f43)


**WHY HISTWEIGHT?**
The MATLAB function histcounts allows to bin scattered data points into quantiles and counts the number of points falling into each quantile.
HISTWEIGHT improves histcounts by considering intensity values for each data point that are spread across the neighbouring quantiles. 
Note that in this way energy conservation is respected as the total sum of the intensity values associated to each point is equal to the total sum of the intensity values associated to each bin.
This does not happen with MATLAB histcounts.
THe following figures show a comparison of HISTWEIGHT against MATLAB histcounts using three different methods for a set of 10 randomly distributed points:
- _invsquared_: inverse squared distance of each vertex with respect to the center of the bin
  ![invsquared](https://github.com/andrepiz/histweight/assets/75851004/f408528c-e8a7-4adb-ab95-c268d866234c)
- _diff_: 1 minus distance of each vertex with respect to the center of the bin, normalized over maximum distance
  ![diff](https://github.com/andrepiz/histweight/assets/75851004/1c50ab79-f9cd-4dfd-bfcb-0f64e9be35aa)
- _area_: fractions of a square box centered into the point that fall into each neighbouring bins
  ![area](https://github.com/andrepiz/histweight/assets/75851004/1209e072-9572-4dc8-9896-73cb8e977fe5)
