from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np

import Paths.Paths as paths
Path = paths.filepaths()


def plot_dendro(fitsfile, table_truncated, table_main,vmin=0,vmax=1):
    

    fitsdata_b6 = fits.open(fitsfile)
    image_b6 = fitsdata_b6[0].data
    if len(image_b6.shape)>2:
        image_b6 = fitsdata_b6[0].data[0][0]
    wcs = WCS(fitsdata_b6[0].header,naxis=2)


    fig = plt.figure(figsize=(30, 30))
    ax = fig.add_subplot(projection=wcs)
    ax.imshow(image_b6, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    ax.scatter(table_main['peak_x'], table_main['peak_y'], s=5, color='orange', marker='o', label='Main Sources')
    ax.scatter(table_truncated['peak_x'], table_truncated['peak_y'], s=5, color='cyan', marker='x',)
    plt.savefig(f'pngs/{region}_{band}_final_dendro.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


for region in ['W51-E', 'W51-IRS2']:
    for band in ['B3', 'B6']:
        if region == 'W51-E':
            if band == 'B3':
                fitsfile = Path.w51e_b3_tt0
             
                vmin = -0.00010907209521789237
                vmax = 0.0009166079520288469
            elif band == 'B6':
                fitsfile = Path.w51e_b6_cont
              
                vmin = -0.00031168037547342546
                vmax = 0.002825007797582483
        elif region == 'W51-IRS2':
            if band == 'B3':
                fitsfile = Path.w51n_b3_tt0
              
                vmin =  -0.00010907209521789237
                vmax = 0.0009166079520288469
            elif band == 'B6':
                fitsfile = Path.w51n_b6_cont
               
                vmin = -0.00031168037547342546
                vmax = 0.002825007797582483
        
        table_truncated = Table.read(f'tables/{region}_{band}_truncated_dendro.fits')
        table_main = Table.read(f'tables/{region}_{band}_initial_dendro.fits')

        print(f'out of {len(table_main)}, {len(table_truncated)} sources are selected for {region} {band} band dendrogram')
        plot_dendro(fitsfile, table_truncated, table_main, vmin=vmin, vmax=vmax)