# Hypsograms

A program to facilitate hypsometric analysis and creation of depth-area metrics for ocean
basins, for use in geochemical and climate modelling.

Currently under development by Morgan Wolfe and Uli Wortmann.
Stay Tuned!

********
## Usage
### Defining a region
```python
my_region = hypsograms.region(name="New Region",bounds = [0,30,0,30],res='10m')
# retrieve elevation data from gmt
my_region.data_from_gmt()
# set elevation limits & compute pixel areas
my_region.clip_z(-9000,0)
my_region.compute_area()
```
### Creating a hypsogram from a region
```python
my_hypsogram = my_region.build_hyps(dz=10)
# get planar area of ocean at 100, 1000, 3000 meters deep 
areas = my_hypsogram.area_at_z([-100,-1000,-3000])
# get continental shelf area
shelf = my_hypsogram.area_between([(-250,0)])
# plot a hypsogram
my_hypsogram.plot()
```
