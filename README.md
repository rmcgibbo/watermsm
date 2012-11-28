watermsm
========

To do
-----

Griding

1. Get the translation to center the CNT atoms
2. Get the rotation to align the CNT atoms' moments of inertia with the
   cartesian axes.
3. Apply this translation and rotation to all of the atoms (water + CNT)
4. Get the number of water oxygens in each grid cell.
   - If the grid width is 1 unit, then griding out the points just amounts to
   rounding the xyz coordinates of each atom to the nearest unit.

   - The next step is one-dimensionalizing the rounded coordinates. i.e. assign
   each of the grid points in R^3 a number in R, probably just
   `(x_i - min(x)) * width(y) * width(z) + (y_i - min(y)) * width(z) + (z_i - min(z)`

   - Then the counting is just histogramming these "flattened" coordinates where
   the bins are just all possible values -- like `np.histogram(a, np.arange(np.max(a)))`
   
   - Now `counts[k]` is the counts in the grid cell k, which is at `z_k = k % width(z)`,
   `y_k = (k - z_k) % width(y)`, `x = k - y_k*width(z) - z) % (width(z)*width(y))`

   - I think this is right.
  
Metric
------

The distance metric is then some vector metric on the counts (histogram) vector.
