

from regions import Regions
import numpy as np
from regions import PixCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from regions import CirclePixelRegion
import astropy.units as u
from radio_beam import Beam

import Paths.Paths as paths
Path = paths.filepaths()

def add_count_column(main_table, table1, table2, table3):
    

    ids1 = set(table1['_name'])
    ids2 = set(table2['_name'])
    ids3 = set(table3['_name'])            

    counts = []
    for row in main_table:
        id_val = row['_name']
        count = int(id_val in ids1) + int(id_val in ids2) + int(id_val in ids3)
        counts.append(count)

    
    main_table.add_column(counts, name='id_count')

    return main_table

labels = ['adam', 'nazar', 'taehwa']
for region in ['W51-E', 'W51-IRS2']:
    for band in ['B3', 'B6']:

        init_dendro = Table.read(f'tables/{region}_{band}_initial_dendro.fits')

        tab1 = Table.read(f'tables/{region}_{band}_{labels[0]}_selected.fits')
        tab2 = Table.read(f'tables/{region}_{band}_{labels[1]}_selected.fits')
        tab3 = Table.read(f'tables/{region}_{band}_{labels[2]}_selected.fits')

        count_added_dendro = add_count_column(init_dendro, tab1, tab2, tab3)

        select_ind = np.where(count_added_dendro['id_count']>=2)
        select_ind2 = np.where(count_added_dendro['id_count']==1)
        print(len(select_ind[0]), len(select_ind2[0]))
        
        truncated_dendro = count_added_dendro[select_ind]
        ambiguous_dendro = count_added_dendro[select_ind2]
        truncated_dendro.write(f'tables/{region}_{band}_truncated_dendro.fits', overwrite=True)
        ambiguous_dendro.write(f'tables/{region}_{band}_ambiguous_dendro.fits', overwrite=True)