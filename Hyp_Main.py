"""
A program to provide spline representations of hypsometric curves for ocean
basins, for the purposes of integration with the ESBMTK
"""
import pygmt
import numpy as np
import pandas as pd
import xarray
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate as interp


def tab2grd(fn: str, outfn: str = "mask_grid.nc"):
    """
    Receives a .tab file and converts it to a NetCDF grid to be used by PyGMT.
    Primarily intended for importing external mask files to the program.

    Parameters
    ----------
    fn : str
        Input filename.
    outfn : str, optional
        Output filename. The default is 'mask_grid.nc'.

    Returns
    -------
    focus_region : list
        List containing the coordinates definining the bounding rectangle of
        the imported grid.
    newgrid : TYPE
        The new pygmt grid created from the input file.

    """
    maskxyz = pd.read_table(fn)

    focus_region = [maskxyz.X.min(), maskxyz.X.max(), maskxyz.Y.min(), maskxyz.Y.max()]

    newgrid = pygmt.xyz2grd(
        x=maskxyz.X,
        y=maskxyz.Y,
        z=maskxyz.MASK,
        region=focus_region,
        spacing=2,
        registration="g",
        outgrid="mask_grid.nc",
    )
    return (focus_region, newgrid)


use_mask = "IO_2d"

focus_region, newgrid = tab2grd(f"{use_mask}.tab")

start = -9000  # lowest elevation
stop = 0  # highest elevation
dz = 10  # height interval in meters

fr = np.round(focus_region)


def get_res_codes(do_print=True):
    """
    Prints the list of shortened strings used by PyGMT to represent map and
    data spatial resolutions. These values are used when communicating with
    PyGMT to obtain data at a variety of sizes. The single letters at the end
    of each string represent the unit - 'd' for degrees; 'm' for minutes; 's'
    for seconds. Can optionally return the list of accepted strings as a 'set'
    object, instead of printing them.

    Parameters
    ----------
    do_print : bool, optional
        Defines whether the resolution codes should be printed. If False,
        the strings will be returned as a set instead. The default is True.

    Returns
    -------
    res_codes : set
        The res codes that are returned if do_print is False.

    """
    res_codes = {
        "01d",
        "30m",
        "20m",
        "15m",
        "10m",
        "06m",
        "05m",
        "04m",
        "03m",
        "02m",
        "01m",
        "30s",
        "15s",
        "03s",
        "01s",
    }
    if do_print:
        print(res_codes)
    else:
        return res_codes


display_region = [-180, 180, -90, 90]

fig = pygmt.Figure()


def get_gmt_data(res="01d", reg=None):
    """
    Uses PyGMT to obtain elevation data at a given resolution, for a given
    region.

    Parameters
    ----------
    res : str, optional
        The string representing the desired resolution of the data. String must
        be formatted in the correct PyGMT style. For a list of accepted values,
        use function get_res_codes(). The default is '01d'.
    reg : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    rel_grid : TYPE
        DESCRIPTION.
    """
    rel_grid = pygmt.datasets.load_earth_relief(resolution=res, region=reg)
    return rel_grid


rel_grid = get_gmt_data(res="01m", reg=fr)

xyz = pygmt.grd2xyz(grid=rel_grid)

xyz_sel = pygmt.select(
    data=xyz, gridmask="mask_grid.nc", z_subregion=str(start) + "/" + str(stop)
)

min_z = xyz_sel.elevation.min()
r_min_z = int(np.round(min_z, -2))

max_z = xyz_sel.elevation.max()
r_max_z = int(np.round(max_z, -2))

xyz_grid = pygmt.xyz2grd(
    data=xyz_sel,
    # x=xyz['lon'],
    # y=xyz['lat'],
    # z=xyz['elevation'],
    spacing="01m",
    region=fr,
    registration="p",
)

fig.coast(
    region=display_region,
    projection="Cyl_stere/12c",
    land="darkgrey",
    water="white",
    shorelines="1/0.5p",
    frame="ag",
)
fig.grdimage(
    grid=xyz_grid,
    projection="Cyl_stere/12c",
    region=display_region,
    nan_transparent=True,
)

fig.coast(
    region=display_region,
    projection="Cyl_stere/12c",
    land=None,
    water=None,
    shorelines="1/0.5p",
    frame="ag",
)
fig.plot([(fr[0], fr[3]), (fr[1], fr[3])], straight_line=True, pen="2p,red")
fig.plot([(fr[0], fr[3]), (fr[0], fr[2])], straight_line=True, pen="2p,red")
fig.plot([(fr[0], fr[2]), (fr[1], fr[2])], straight_line=True, pen="2p,red")
fig.plot([(fr[1], fr[3]), (fr[1], fr[2])], straight_line=True, pen="2p,red")


fig.show()


def geo_radius(lat=None):
    """
    Returns the approximate radius of the earth for a given latitude.

    Parameters
    ----------
    lat: TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    R: TYPE
        DESCRIPTION.

    """
    r_lat = np.radians(lat)
    # R = 6378000 * np.sqrt((1-(2*(e**2)-(e**4))*np.sin(r_lat)**2) /
    #                       (1-((e**2)*np.sin(r_lat)**2)))
    a: int = 6378137  # equatorial radius
    b: int = 6356752  # polar radius

    R: float = (
        ((a ** 2 * np.cos(r_lat)) ** 2 + (b ** 2 * np.sin(r_lat)) ** 2)
        / ((a * np.cos(r_lat)) ** 2 + (b * np.sin(r_lat)) ** 2)
    ) ** (1 / 2)

    return R


def grid_area(lat: float, lon: float, size: float) -> float:
    """Calculate the area of a rectangular area of size = 1 deg at a
    given lat long position

    """

    r: float = geo_radius(lat) / 1000
    dy: float = (size * r * np.pi) / 180
    dx: float = size / 180 * np.pi * r * np.cos(np.radians(lat))

    a: float = abs(dx * dy)
    return a


res = 1  # id resolution in minutes
r = f"{res}m"
grid = xyz_grid

dx = res / 60  # convert resolution in degrees
steps = int(round(abs(fr[2] - fr[3]) / dx))


# initialize


def init_hyp():
    lat = np.linspace(fr[2] + dx, fr[3] - dx, steps)
    weight = np.linspace(fr[2] + dx, fr[3] - dx, steps)
    a = np.linspace(fr[2] + dx, fr[3] - dx, steps)
    count = np.arange(start, stop, dz, dtype=float)
    elevation = np.arange(start, stop, dz)
    return (lat, weight, a, count, elevation)


lat, weight, a, count, elevation = init_hyp()

""" calc weighing factor as function of latitude. If we ignore that earth
is not a perfect sphere, this would simply be the cosinus of the latidude
"""
# for i, l in enumerate(lat):
#     weight[i] = grid_area(l, 0, dx)

weight = grid_area(lat, 0, dx)

""" Loop over elevation intervals, and count occurances.
Adjust the count for each latidude to reflect the changing
data density
"""
for i, e in enumerate(elevation):
    # this will return a latitudinal transect
    a = np.sum(np.logical_and(grid > e, grid < e + dz), axis=1)
    count[i] = np.sum(a * weight)
# normalize cout per transect to total number of counts
count = count / np.sum(count)
# get the cumulative count distribution
cum = np.cumsum(count)
df = pd.DataFrame({"Elevation": elevation, "CumSum": cum, "Count": count})
df.to_csv(f"{use_mask}_{start}_{stop}_Hypsometric_Curve_{r}.csv", float_format="%8.32f")

# check plot
fig1, ax1 = plt.subplots()  #
ax1.plot(cum, elevation)
ax1.set_ylabel("Elevation [m]")
ax1.set_xlabel(f"Area  [cumulative %]")
ax1.invert_xaxis()
plt.show(block=False)

splin = interp.splrep(elevation, cum)

splin_t, splin_c, splin_k = splin
