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
