# -*- coding: utf-8 -*-
"""
Tabulates provided data for specified ocean basins and returns the result as
a Pandas DataFrame.
"""
import pygmt
import numpy as np
import pandas as pd
import xarray
import matplotlib.pyplot as plt

cats = {'Trench': (-9000, -6000),
        'Basin': (-6000, -3000),
        'Slope': (-3000, -250),
        'Shelf': (-250, 0)}


def tabulate_data(categories: dict = cats):
    pass


elev_range = '-9000_0'
pacific = pd.read_csv(f"PO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
atlantic = pd.read_csv(f"AO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
indian = pd.read_csv(f"IO_2d_{elev_range}_Hypsometric_Curve_1m.csv")
southern = pd.read_csv(f"SO_2d_{elev_range}_Hypsometric_Curve_1m.csv")

cats = {'Trench': (-9000, -6000),
        'Basin': (-6000, -3000),
        'Slope': (-3000, -250),
        'Shelf': (-250, 0)}

p_tot = np.empty(4)
a_tot = np.empty(4)
i_tot = np.empty(4)
s_tot = np.empty(4)

files = ['PO', 'AO', 'IO', 'SO']


cat_names = ['Trench', 'Basin', 'Slope', 'Shelf']

for i, e in enumerate(cat_names):
    p_tot[i] = np.sum(pacific.Count[cats[e][0] <= pacific.Elevation]
                      [pacific.Elevation <= cats[e][1]])
    a_tot[i] = np.sum(atlantic.Count[cats[e][0] <= atlantic.Elevation]
                      [atlantic.Elevation <= cats[e][1]])
    i_tot[i] = np.sum(indian.Count[cats[e][0] <= indian.Elevation]
                      [indian.Elevation <= cats[e][1]])
    s_tot[i] = np.sum(southern.Count[cats[e][0] <= southern.Elevation]
                      [southern.Elevation <= cats[e][1]])

tab_df = pd.DataFrame({'Pacific': p_tot, 'Atlantic': a_tot,
                      'Indian': i_tot, 'Southern': s_tot}, index=cat_names)
print(tab_df)
