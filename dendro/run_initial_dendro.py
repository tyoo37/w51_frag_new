from dendrocat import RadioSource
from astropy.io import fits
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import sys
import importlib
from astropy.wcs import WCS
from regions import Regions, PixCoord
from astropy import stats
from itertools import chain
import dendrocat
import regions

import Paths.Paths as paths
Path = paths.filepaths()


def remove_dendro(source_object,custom_dendro,pbdat,remove_regs, thres=0.21,):
    
    structure = np.array(custom_dendro.leaves)
    table = source_object.catalog
    
    
  
    xarr = []; yarr=[] 
    peakvalarr = []

    
    for s in structure:
        
        xind = s.indices(subtree=True)[0]
        yind = s.indices(subtree=True)[1]

        value = s.values(subtree=True)
        xpos = yind[np.argmax(value)]
        ypos = xind[np.argmax(value)]
        
        xarr.append(xpos) ; yarr.append(ypos) ; peakvalarr.append(np.max(value)) 
    table.add_column(xarr, name='peak_x')
    table.add_column(yarr, name='peak_y')
    table.add_column(peakvalarr, name='peak_value')

    remove_idx =[]
    if remove_regs is not None:
        for i in range(len(table['peak_x'])):
            peakx = table['peak_x'][i]
            peaky = table['peak_y'][i]
            if pbdat is not None:
                if pbdat[peakx.astype(int), peaky.astype(int)] < thres:
                    remove_idx.append(i)
                    continue
            if remove_regs is not None:
                for reg in remove_regs:
                    if reg.contains(PixCoord(peakx,peaky)):
                        remove_idx.append(i)
                        break
    table.remove_rows(remove_idx)

    return table


    
   

def main():
    """
    
    """
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-r", "--region", dest="region",
                    default='W51E',
                    help="region", metavar="region")
    parser.add_option("--band", dest="band",
                    default='B3',
                    help="band", metavar="band")
    
    (options, args) = parser.parse_args()

    region = options.region
    band = options.band
    print(f'run_initial_dendro.py: Running initial dendrogram for {region} {band}')
    if region == 'W51-E':
        if band == 'B3':
            fitsfile = Path.w51e_b3_tt0
            pbfile = None
            noiseregion = Path.w51e_b3_noise_region
            #regfile_to_remove = Path.w51e_b3_remove_reg
            regfile_to_remove = None
            vmin = -0.00010907209521789237
            vmax = 0.0009166079520288469
        elif band == 'B6':
            fitsfile = Path.w51e_b6_tt0
            pbfile = None
            noiseregion = Path.w51e_b6_noise_region
            regfile_to_remove = None
            vmin = -0.00031168037547342546
            vmax = 0.002825007797582483
    elif region == 'W51-IRS2':
        if band == 'B3':
            fitsfile = Path.w51n_b3_tt0
            pbfile = None
            noiseregion = Path.w51n_b3_noise_region
            #regfile_to_remove = Path.w51n_b3_remove_reg
            regfile_to_remove = None
            vmin =  -0.00010907209521789237
            vmax = 0.0009166079520288469
        elif band == 'B6':
            fitsfile = Path.w51n_b6_tt0
            pbfile = None
            noiseregion = Path.w51n_b6_noise_region
            regfile_to_remove = None
            vmin = -0.00031168037547342546
            vmax = 0.002825007797582483
    fitsdata_b6 = fits.open(fitsfile)
    noiseregion_b6 = Regions.read(noiseregion,format='ds9')
    image_b6 = fitsdata_b6[0].data
    if len(image_b6.shape)>2:
        image_b6 = fitsdata_b6[0].data[0][0]
    wcs = WCS(fitsdata_b6[0].header,naxis=2)

    container = []
    for reg in noiseregion_b6:
        pix_reg = reg.to_pixel(wcs)
        noisemask = pix_reg.to_mask()
        noiseim_b6 = noisemask.cutout(image_b6)
        
        container.append(noiseim_b6.flatten())

    noiseim_b6 = np.concatenate(container)
    b6_std = stats.mad_std(noiseim_b6,ignore_nan=True)

    min_value_factor = 3
    min_delta_factor = 1.5
    min_npix = 15

    source_object_b6 = RadioSource(fitsdata_b6,)


    custom_dendro_b6 = source_object_b6.to_dendrogram(min_value=min_value_factor*b6_std, 
                                                    min_delta=min_delta_factor*b6_std, 
                                                    min_npix=min_npix)
    table_b6 = source_object_b6.to_catalog()



    if regfile_to_remove is not None:
        remove_regs = Regions.read(regfile_to_remove, format='ds9')
    else:
        remove_regs = None
    if pbfile is not None:
        pb_fits = fits.open(pbfile)
        pbdata = pb_fits[0].data
        if len(pbdata.shape)>2:
            pbdata = pb_fits[0].data[0][0]
    else:
        pbdata = None
    table_b6 = remove_dendro(source_object_b6, custom_dendro_b6, pbdata, remove_regs, thres=0.21)
    table_b6.write(f'tables/{region}_{band}_initial_dendro.fits', format='fits', overwrite=True)

    fig = plt.figure(figsize=(30, 30))
    ax = fig.add_subplot(projection=wcs)
    ax.imshow(image_b6, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    ax.scatter(table_b6['peak_x'], table_b6['peak_y'], s=100, color='red', marker='x',)
    plt.savefig(f'pngs/{region}_{band}_initial_dendro.png', dpi=300, bbox_inches='tight')
    plt.close(fig)



if __name__ == "__main__":
    main()