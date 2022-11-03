# -*- coding: utf-8 -*-
"""
Hypsogram.

"""


class hypsogram:
    def __init__(self, region, htype):
        # TODO add attributes for a hypsogram instance (per Uli's recs)
        self.region = region  # Atlantic, Pacific, Indian, Southern, etc...
        self.htype = htype  # (lookup, spline)
        # zmax, zmin
        pass

    def area(z):
        # TODO return area(s) at given elevation(s)
        pass

    def area_between(zs):
        # TODO return area(s) between given elevation intervals
        pass

    def volume():
        # TODO return volume between given elevation interval(s)
        pass

    def make_spline():
        # TODO interpolate spline fn parameters from tabled data
        pass

    def make_table():
        # TODO generate a table of fixed values from spline parameters
        # unsure when this would be used, but low effort to implement
        # & keeps program balanced
        pass


class region:
    def __init__(self):
        # TODO add attributes to define existing oceans calculated
        # name
        # filepath (if applic)
        # max elevation range defined
        # TODO define existing regions as part of library setup
        pass


# ??? pass display of hypsometric curves to mpl? or leave as exercise for user?
