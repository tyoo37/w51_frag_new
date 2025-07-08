from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from regions import CircleSkyRegion
from regions.core import PixCoord  # Not used here, but good to know
import Paths.Paths as paths
import numpy as np
Path = paths.filepaths()

def cat_to_crtf(catfile, output_file, radius=0.04 * u.arcsec, color='cyan', width=2, no_use_ra=False):
    # Step 1: Load your table
    print(catfile)
    tbl = Table.read(catfile)  # Or .fits, .csv, etc.

    # Step 2: Create SkyCoord from the table
    if not no_use_ra:
        coords = SkyCoord(ra=tbl['ra'] * u.deg, dec=tbl['dec'] * u.deg, frame='icrs')
    else:
        ra = []
        dec = []
        for i in range(len(tbl['b3_xsky'])):
            if tbl['b6_xsky'][i] < 0:
                ra.append(tbl['b3_xsky'][i]) 
                dec.append(tbl['b3_ysky'][i])
            else:
                ra.append(tbl['b6_xsky'][i]) 
                dec.append(tbl['b6_ysky'][i])
        coords = SkyCoord(ra=np.array(ra) * u.deg, dec=np.array(dec) * u.deg, frame='icrs')

    # Step 3: Create CRTF region strings
    region_strings = []
    for i, coord in enumerate(coords):
        # Format: circle [[ra, dec], radius]
        region_name=str(i)
        region_string = f"circle [[{coord.ra.deg}deg, {coord.dec.deg}deg], {radius.to(u.arcsec).value}arcsec] coord=ICRS color={color} width={width} label={region_name}"
        region_strings.append(region_string)

    # Step 4: Write regions to file
    with open(output_file, 'w') as f:
        f.write("#CRTFv0\n")
        f.write("\n".join(region_strings))

cat_to_crtf(Path.w51e_dendro_matched_catalog, 'tables/w51e_dendro_matched_catalog_old.crtf', color='blue', no_use_ra=True)
cat_to_crtf(Path.w51n_dendro_matched_catalog, 'tables/w51n_dendro_matched_catalog_old.crtf', color='green', no_use_ra=True)
cat_to_crtf('tables/dendro_w51e_matched.fits', 'tables/w51e_dendro_matched_catalog_new.crtf', color='red', radius=0.07 * u.arcsec)
cat_to_crtf('tables/dendro_w51n_matched.fits', 'tables/w51n_dendro_matched_catalog_new.crtf', color='orange', radius=0.07 * u.arcsec)
cat_to_crtf('tables/dendro_w51e_matched_duplicated_removed.fits', 'tables/w51e_dendro_matched_catalog_new_duplicated_removed.crtf', color='magenta', radius=0.1 * u.arcsec)