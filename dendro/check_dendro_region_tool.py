
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


def indexing_regions(regfile):
    reg_list = []
    regs = Regions.read(regfile, format='ds9')
    with open(regfile, 'r') as f:
        for idx, line in enumerate(f):
            if 'width' in line:
                reg_list.append(line)

    b3b6_regs = []
    b3only_regs = []
    for i, reg in enumerate(reg_list):
        if 'color=#cyan' in reg:
            b3b6_regs.append(i)
        elif 'color=#red' in reg:
            b3only_regs.append(i)
    return b3b6_regs, b3only_regs



def compare_with_initial_dendro(init_dendro, regs, wcs, beam=None, ispointregion=False, label='adam'):
    pixcoords = PixCoord(init_dendro['peak_x'], init_dendro['peak_y'])

    selected_regs = []
    unselected_regs = []
    selected_peaks_x = []
    selected_peaks_y = []
    idarr = np.arange(len(init_dendro['peak_x']))
    selected_idarr = []

    pixel_scale = wcs.proj_plane_pixel_scales()[0]
    copied_dendro = init_dendro.copy()
    if not f'{label}_selected' in init_dendro.colnames:
        copied_dendro.add_column(np.zeros(len(init_dendro)), name=f'{label}_selected')

    for reg in regs:
        reg_pix = reg.to_pixel(wcs)
        if ispointregion:
            reg_pix = CirclePixelRegion(reg_pix.center, radius=(beam.major/pixel_scale).to(u.pix/u.pix).value)

        if any(reg_pix.contains(pixcoords)):
           # print('matched with dendro')
            selected_regs.append(reg)

            iscontain = reg_pix.contains(pixcoords)
            selected_peak_x = init_dendro['peak_x'][iscontain]
            selected_peak_y = init_dendro['peak_y'][iscontain]
            selected_peak_value = init_dendro['peak_value'][iscontain]
            selected_id = idarr[iscontain]
            if len(selected_peak_x) > 1:
                maxind = np.argmax(selected_peak_value)
                selected_peaks_x.append(selected_peak_x[maxind].item())
                selected_peaks_y.append(selected_peak_y[maxind].item())
                selected_idarr.append(selected_id[maxind].item())
                
                #copied_dendro[f'{label}_selected'][iscontain][maxind] = 1

            elif len(selected_peak_x) == 1:
                selected_peaks_x.append(selected_peak_x.item())
                selected_peaks_y.append(selected_peak_y.item())
                selected_idarr.append(selected_id.item())
                #copied_dendro[f'{label}_selected'][iscontain] = 1



        else:
#            print('not matched with dendro')
            unselected_regs.append(reg)
    print(f'Number of selected regions: {len(selected_regs)}')
    print(f'Number of unselected regions: {len(unselected_regs)}')
    print(f'Number of total regions: {len(regs)}')
    print(f'Number of selected initial dendrogram peaks: {len(selected_peaks_x)}')
    print(f'Number of initial dendrogram peaks: {len(init_dendro)}')
    print(selected_peaks_x)
    selected_dendro = init_dendro[selected_idarr]
    return selected_regs, unselected_regs, selected_peaks_x, selected_peaks_y, selected_dendro



def plot_selected_regions(ax, wcs , selected_regs, unselected_regs, selected_peaks_x, selected_peaks_y, x, y, ispointregion=False, beam=None):
    from regions import Regions
    from astropy.visualization.wcsaxes import WCSAxes

    pixel_scale = wcs.proj_plane_pixel_scales()[0]
    

    if ispointregion:
        for reg in selected_regs:
            reg_pix = reg.to_pixel(wcs)
            reg_pix_circle = CirclePixelRegion(reg_pix.center, radius=(beam.major/pixel_scale).to(u.pix/u.pix).value)
            reg_pix_circle.plot(ax=ax,  color='cyan', lw=1)
    
        for reg in unselected_regs:
            reg_pix = reg.to_pixel(wcs)
            reg_pix_circle = CirclePixelRegion(reg_pix.center, radius=(beam.major/pixel_scale).to(u.pix/u.pix).value)
            reg_pix_circle.plot(ax=ax,  color='red', lw=1)
    else:    
        for reg in selected_regs:
            reg_pix = reg.to_pixel(wcs)
            reg_pix.plot(ax=ax, facecolor='none', edgecolor='cyan', lw=1)
        
        for reg in unselected_regs:
            reg_pix = reg.to_pixel(wcs)
            reg_pix.plot(ax=ax, facecolor='none', edgecolor='red', lw=1)
    ax.scatter(x, y, marker='x', color='orange')
    ax.scatter(selected_peaks_x, selected_peaks_y, color='blue', marker='x', label='Peaks')



def plot_main(region, band, fitsfile, init_dendrofile, regs, label='adam', ispointregion=False):
    

    fits_data = fits.open(fitsfile)
    hdr = fits.getheader(fitsfile)
    init_dendro = Table.read(init_dendrofile)
    wcs = WCS(hdr, naxis=2)
    if ispointregion:
        beam = Beam.from_fits_header(hdr)
    else:
        beam = None
    selected_regs, unselected_regs, selected_peaks_x, selected_peaks_y, copied_dendro = compare_with_initial_dendro(init_dendro, regs, wcs, beam=beam, ispointregion=ispointregion, label=label)
    print(selected_peaks_x, selected_peaks_y)
    fig = plt.figure(figsize=(30, 30))
    ax = fig.add_subplot(111, projection=WCS(hdr, naxis=2))
    ax.imshow(fits_data[0].data[0][0], origin='lower', cmap='gray', vmin=0, vmax=1.5e-3)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    plot_selected_regions(ax, wcs, selected_regs, unselected_regs, selected_peaks_x, selected_peaks_y, init_dendro['peak_x'], init_dendro['peak_y'], ispointregion=ispointregion, beam=beam)
    plt.savefig(f'pngs/{region}_{band}_{label}_selected_regions.png', bbox_inches='tight', dpi=300)

    return copied_dendro

def save_dendro(selected_dendro, region, band, label='adam'):
    
    
    selected_dendro.write(f'tables/{region}_{band}_{label}_selected.fits', overwrite=True)
    print(f'Saved selected dendrogram to tables/{region}_{band}_{label}_selected.fits')