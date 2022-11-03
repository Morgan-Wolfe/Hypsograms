# -*- coding: utf-8 -*-
"""
Plots hypsographic curves from provided datasets.
"""
import pandas as pd
import matplotlib.pyplot as plt


elev_range = "-9000_0"

fig, axmain = plt.subplots(1, 2, sharey=True)
td = pd.read_csv(f"PO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
axmain[0].plot(td.Count, td.Elevation, label="Pacific Ocean")
axmain[1].plot(td.CumSum, td.Elevation, label="Pacific Ocean")

td = pd.read_csv(f"Globe_2d_{elev_range}_Hypsometric_Curve_1m.csv")
axmain[0].plot(td.Count, td.Elevation, label="Global Curve")
axmain[1].plot(td.CumSum, td.Elevation, label="Global Curve")

td = pd.read_csv(f"IO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
axmain[0].plot(td.Count, td.Elevation, label="Indian Ocean")
axmain[1].plot(td.CumSum, td.Elevation, label="Indian Ocean")

td = pd.read_csv(f"AO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
axmain[0].plot(td.Count, td.Elevation, label="Atlantic Ocean")
axmain[1].plot(td.CumSum, td.Elevation, label="Atlantic Ocean")

td = pd.read_csv(f"SO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
axmain[0].plot(td.Count, td.Elevation, label="Southern Ocean")
axmain[1].plot(td.CumSum, td.Elevation, label="Southern Ocean")

axmain[0].legend()
axmain[0].grid()
axmain[0].set_title("% Area per 10m")
axmain[0].hlines([-9000, -6000, -3000, -250], 0, 0.01, color="k")
axmain[0].set_ylabel("Elevation (m)")
axmain[0].set_xlabel("Area (%)")

axmain[1].set_xlabel("Cumulative Area (%)")
axmain[1].invert_xaxis()
axmain[1].legend()
axmain[1].set_title("Hyspographic Curve")
axmain[1].grid()
axmain[1].hlines([-9000, -6000, -3000, -250], 0, 1, "k")

fig.subplots_adjust(wspace=0.0)
fig.suptitle("Ocean Basin Hypsographies and Area Distributions by Elevation")
fig.show()
