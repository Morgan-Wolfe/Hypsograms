import numpy as np
import pandas as pd
import xarray
import scipy
import scipy.interpolate as interp
import pygmt
import matplotlib.pyplot as plt


def grid_area(lat: float, lon: float, size: float = 1) -> float:
    """
    Calculate the area of a rectangular area of size = 1 deg at a
    given lat long position

    Parameters
    ----------
    lat : float
        DESCRIPTION.
    lon : float
        DESCRIPTION.
    size : float, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    float
        Represented area in square kilometres.

    """
    r_lat = np.radians(lat)
    r_lon = np.radians(
        lon
    )  # alternate calc'n utilizing lat & lon anticipated in future
    # R = 6378000 * np.sqrt((1-(2*(e**2)-(e**4))*np.sin(r_lat)**2) /
    #                       (1-((e**2)*np.sin(r_lat)**2)))

    a: int = 6378137  # equatorial radius
    b: int = 6356752  # polar radius

    R: float = (
        ((a**2 * np.cos(r_lat)) ** 2 + (b**2 * np.sin(r_lat)) ** 2)
        / ((a * np.cos(r_lat)) ** 2 + (b * np.sin(r_lat)) ** 2)
    ) ** (
        1 / 2
    )  # calculate earth radius at given latitude

    r: float = R / 1000  # convert to km
    dy: float = (size * r * np.pi) / 180  # calc long'l side length
    dx: float = (
        size / 180 * np.pi * r * np.cos(np.radians(lat))
    )  # calc lat side

    a: float = abs(dx * dy)  # area
    return a


def test_random_region(res="01d"):
    """For testing/debug.

    Creates a region with random bounds and a resolution of 1 degree as a
    sample for function and capabilities tests.

    Parameters
    ----------
    res : TYPE, optional
        DESCRIPTION. The default is "01d".

    Returns
    -------
    reg : TYPE
        DESCRIPTION.

    """
    x = np.random.randint(-180, 180, 2)
    x.sort()
    y = np.random.randint(-90, 90, 2)
    y.sort()
    bounds = np.append(x, y)
    reg = region(bounds=bounds, res=res)
    reg.data_from_gmt(bounds=bounds, res=res)
    return reg


gmt_rescodes = {
    "01d": 1,
    "30m": 30 / 60,
    "20m": 20 / 60,
    "15m": 15 / 60,
    "10m": 10 / 60,
    "06m": 6 / 60,
    "05m": 5 / 60,
    "04m": 4 / 60,
    "03m": 3 / 60,
    "02m": 2 / 60,
    "01m": 1 / 60,
    "30s": 0.5 / 60,
    "15s": 0.25 / 60,
    "03s": 0.05 / 60,
    "01s": 1 / 60 / 60,
}


class region:
    """
    Main class defining a geographic region as a rectangular bounding box.
    corners defined by region.bounds
    """

    def __init__(
        self,
        name="New Region",  # name - future plans
        bounds=[0, 0, 0, 0],
        data=None,
        res="01d",
        parent=None,
    ):
        self.bounds = np.round(bounds)
        self.name = name
        self.data = data
        self.area = None
        self.mdat = None
        self.interval = gmt_rescodes[res]
        self.parent = parent
        self.res = res
        # TODO add attributes to define existing oceans calculated
        # TODO define existing regions as part of library setup
        pass

    def build_hyps(self, dx=None, dz=10, name=None):
        """
        Constructs a hypsogram from the parent region.

        Parameters
        ----------
        dx : int, optional
            DESCRIPTION. The default is None.
        dz : TYPE, optional
            DESCRIPTION. The default is 10m.
        name : TYPE, optional
            DESCRIPTION. Defaults to the same name as region.

        Returns
        -------
        Hypsogram
            DESCRIPTION.

        """
        if name == None:
            name = self.name
        if dx == None:
            dx = self.interval
        grid = self.data.elevation

        elevation = np.arange(
            int(self.data.elevation.min()), int(self.data.elevation.max()), dz
        )
        count = np.arange(
            int(self.data.elevation.min()), int(self.data.elevation.max()), dz
        )
        area = np.arange(
            int(self.data.elevation.min()), int(self.data.elevation.max()), dz
        )
        weight = grid_area(self.data.lat, self.data.lon, dx)
        for i, e in enumerate(elevation):
            # this will return a latitudinal transect
            a = np.sum(np.logical_and(grid > e, grid < e + dz), axis=1)
            area[i] = np.sum(a * weight)
        # normalize cout per transect to total number of counts
        count = area / np.sum(area)
        cum = np.cumsum(count)
        ca = np.cumsum(area)
        # dftest = pd.DataFrame(
        #     {
        #         "Elevation": elevation,
        #         "CumSum": cum,
        #         "Count": count,
        #         "Area": area,
        #         "CumArea": ca,
        #     }
        # )
        # vals = ["CumSum", "Count", "Area"]
        # datest = xarray.DataArray(
        #     [cum, count, area],
        #     coords=[elevation, vals],
        #     dims=["Elevation", "Data"],
        # )

        ds = xarray.Dataset(
            {
                "CumSum": (["Elevation"], cum),
                "Count": (["Elevation"], count),
                "Area": (["Elevation"], area),
                "CumArea": (["Elevation"], ca),
            },
            coords={"Elevation": (["Elevation"], elevation)},
        )

        return Hypsogram(region=self, htype=1, data=ds, zstep=dz)

    def data_from_gmt(self, bounds=None, res=None):
        """retrieves data from GMT at the provided location and resolution"""
        if res is None:
            res = self.res
        elif res in gmt_rescodes.keys():
            self.res = res

        self.interval = gmt_rescodes[self.res]

        if bounds is None:
            bounds = self.bounds
        rel_grid = pygmt.datasets.load_earth_relief(
            resolution=res, region=bounds
        )
        self.data = rel_grid

    def compute_area(self):
        """
        Calculates approximate planar area represented by each datapoint within
        the region. Current method assumes no variation in Earth radius by lon.
        The area represented by each 'cell' is then added to the region data as
        a new dimension.

        Returns
        -------
        None.

        """
        if self.area is not None:
            self.area = None
        area_grid = grid_area(self.data.lat, self.data.lon, self.interval)
        area_grid = area_grid * self.data.notnull()
        area_grid.name = "area"
        self.data = xarray.merge([self.data, area_grid])
        self.area = float(area_grid.sum())

    def plot_map(self, proj="aitoff"):
        """
        Displays a map illustrating the location and density of the region's
        non-null datapoints.

        Parameters
        ----------
        proj : String, optional
            Represents the map projection to be used. The default is "aitoff".

        Returns
        -------
        None.

        """
        datastack = self.data.stack(coords=("lat", "lon"))
        nonulls = datastack.where(datastack.isnull() == False, drop=True)
        plt.figure()
        plt.subplot(projection=proj)
        plt.scatter(
            (nonulls.lon / 180) * np.pi, (nonulls.lat / 90) * (np.pi / 2), s=1
        )
        plt.grid(True)

    def clip_z(self, zmin=None, zmax=0):
        """
        Sets elevations within the region that fall outside the given range to
        null values. Often used to clip land (positive) elevations from marine
        data.

        Parameters
        ----------
        zmin : Float, optional
            DESCRIPTION. Defaults to the region's lowest elevation.
        zmax : Float, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        """
        if zmin is None:
            zmin = self.data.min()
        self.data = self.data.where(
            np.logical_and(zmin <= self.data, self.data < zmax)
        )


class Hypsogram:
    def __init__(
        self, region, htype, data=None, atype=1, name="Hypsogram", zstep=None
    ):
        self.region = region  # Atlantic, Pacific, Indian, Southern, etc....
        self.htype = htype  # (lookup, spline)
        self.data = data
        self.zmin = data.Elevation.min()
        self.zmax = data.Elevation.max()
        self.atype = atype
        self.area = region.area
        self.mdata = {"reg_meta": region.mdat}
        self.name = name
        self.zstep = zstep
        pass

    def area_at_z(self, z, isPlanar: bool = True):
        """
        Returns the areas for the given list of elevations.

        Parameters
        ----------
        z : List
            DESCRIPTION.
        isPlanar : bool, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        areas : List
            DESCRIPTION.

        """
        areas = []
        for i in z:
            ap = self.data.Area.sel(Elevation=i, method="nearest")
            areas.append(ap)
        return areas

    def area_between(self, zs, isPlanar: bool = True):
        """
        Returns the combined area between the elevations given.

        Parameters
        ----------
        zs : List of tuples
            DESCRIPTION.
        isPlanar : bool, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        areas : List
            DESCRIPTION.

        """
        areas = []
        for i in zs:
            areas.append(
                self.data.Area.where(
                    np.logical_and(
                        self.data.Elevation >= i[0],
                        self.data.Elevation < i[1],
                    )
                ).sum()
            )
        return areas

    def volume_between(self, z1=None, z2=None):
        """
        Returns the volume between two elevations

        Parameters
        ----------
        z1 : TYPE, optional
            DESCRIPTION. The default is None.
        z2 : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        volume : float
            Volume, in cubic km.

        """
        if z1 == None:
            z1 = self.zmin
        if z2 is None:
            z2 = self.zmax
        h = abs(z1 - z2) / 1000
        area = self.area_between([(z1, z2)])
        volume = h * area
        return volume

    def plot(self, mode="standard"):
        """
        Plots a hyposgram, or elevation-frequency graph, depending on mode.

        Parameters
        ----------
        mode : string, optional
            DESCRIPTION. The default is "standard".

        Returns
        -------
        None.

        """
        import matplotlib.pyplot as plt

        if mode == "standard":
            plt.plot(self.data.CumArea, self.data.Elevation)
        elif mode == "frequency":
            plt.plot(self.data.Area, self.data.Elevation)
