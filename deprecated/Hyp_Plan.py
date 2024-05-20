# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 11:41:40 2022

@author: ibmot
"""

# =============================================================================
# class Region
# #Defines an area where calculations and comparisons will be performed
#
# region.init
#
#
#
# region.area
#
# region.volume
#
# region.bounds
#
# region.subregions
#
# region.newsubregion
#
#
# class Bounds
# #A tuple of coordinates, or a list of tuples of coordinates defining a region
#
# =============================================================================


class Basemap:
    """Defines base elevation data used for hypsograpy. Can be generated
    automatically from GMT or provided by the user."""

    def __init__(self):
        # self.bounds
        # self.resolution
        # self.elevation
        #

        pass

    pass


class Bounds:
    """
    Composed of a tuple of coordinates, or a list of tuples of coordinates,
    defining a geographic area.
    """

    # should be able to test if point within or outside bounds
    # should be able

    def __init__(self):
        pass

    def within(self, bounds):
        "Is self within bounds"
        pass

    pass


class Region:
    """
    The region of interest that is to be calculated upon. Some calcs/methods
    would only need self, others may require additional Regions
    """

    # area defined by a Bounds
    # each Region instance can only have a single Bounds
    # Regions should either contain or refer to the data they relate
    # Regions can have subregions
    #

    def __init__(self, bounds):
        self.bounds = bounds

    pass


class Basin:
    def __init__(self, name):
        self.name = name


# -------- Functions
#       - display hypsographic curve
#       - generate hypsographic splines
#       - display elevation map
#       - summarize
#       - compare
#       - split/clip elevations
#       - stack/merge elevations
#       - calculate contributions
#       -

# ========================
# find name for library             -   DONE
# set up github acct and project    -   DONE
# share w Uli
#
# ===================
# TODO from Uli
# The library should have a classes that returns a hypsometric lookup table instance

# import hypsography as hpy

# A.hyp = hpy(region="Atlantic",
#         h_type='lookuptable', # or spline_approximation)
#         max_z = 200  # relative to sealevel
#         min_z = -8000 # relative to sealevel
#         )
# where the type can be either an integer based lookup table where the index = elevation in meters, or a function that calculates the relevant data from a spline approximation. Min and max z are there to limit the size of the lookup table (no need to carry the upper 8000 meters if one does the marine stuff). All of these should be based on pre-evaluated data (i.e. the spline parameters exist already (see below), so the initialization is quick

# This instance must have the following methods

# A.hyp.area_z(z) = # area at a given elevation
# A.hyp.area_dz(max,min) = # area between two elevations
# A.hyp.volume_dz(max,min) = # volume between two elevations
# Additional useful methods could be:

# generate spline approximation from tabled data
# boot strap hypsometric curve for a given region (this is the least important function, since one would typically only do this, if there is a major correction to the elevation data).
