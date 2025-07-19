from regions import Regions
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
from collections import defaultdict
import numpy as np
import networkx as nx
from regions import CircleSkyRegion
from regions import PixCoord
from astropy.coordinates import search_around_sky
import Paths.Paths as paths
Path = paths.filepaths()

def check_overlap(image_dir1, image_dir2, cat1, cat2):
    cats = [cat1, cat2]
    catnum1 = len(cat1['peak_x_sky'])
    catnum2 = len(cat2['peak_y_sky'])
    isinside1 = np.zeros(catnum1, dtype=bool)
    isinside2 = np.zeros(catnum2, dtype=bool)
    isinsides = [isinside1, isinside2]
    image_dirs = [image_dir1, image_dir2]
    catnums = [catnum1, catnum2]
    for i, image_dir in enumerate(image_dirs):
        fitsdata = fits.open(image_dir)
        wcs = WCS(fitsdata[0].header,naxis=2)
        image = fitsdata[0].data[0][0]
        if len(image.shape)>2:
            image = fitsdata[0].data[0][0]

        cat_comp = cats[i] ; cat = cats[i-1]

        xpeak = cat['peak_x_sky']
        ypeak = cat['peak_y_sky']

        xypeak = np.vstack((xpeak,ypeak)).T

        xypix = wcs.wcs_world2pix(xypeak,0)
        
        for j in range(catnums[i-1]):
            if xypix[j,0] < 0 or xypix[j,1] < 0 or xypix[j,0]>=image.shape[0] or xypix[j,1] >= image.shape[1]:
                isinsides[i-1][j] = False
            
            else:
                isinsides[i-1][j] = np.isfinite(image[int(xypix[j,0]),int(xypix[j,1])])
                
    return isinsides



def get_insignificant_regions(region_a_file, region_b_file, region_c_file, tolerance=1.1e-5*u.deg):
    regions_A = Regions.read(region_a_file)
    regions_B = Regions.read(region_b_file)
    regions_C = Regions.read(region_c_file)

    region_lists = [regions_A, regions_B, regions_C]
    
    def get_coords(region_list):
        return SkyCoord([r.center.icrs for r in region_list])

    coords_lists = [get_coords(regs) for regs in region_lists]

    # === Build overlap graph ===
    G = nx.Graph()

    def tag(i, j):
        return f"{i}_{j}"

    # Find all pairs of overlapping regions between different files
    for i in range(len(coords_lists)):
        for j in range(i+1, len(coords_lists)):
            idx1, idx2, _, _ = search_around_sky(coords_lists[i], coords_lists[j], tolerance)
            for a, b in zip(idx1, idx2):
                G.add_edge(tag(i, a), tag(j, b))

    # === Filter components with ≥2 sources ===
    qualified_components = []
    for component in nx.connected_components(G):
        sources = set(n.split("_")[0] for n in component)
        if len(sources) >= 2:
            qualified_components.append(component)

    # === Merge regions from each qualified component ===
    def get_region(source_idx):
        src, idx = map(int, source_idx.split("_"))
        return region_lists[src][idx]

    combined_regions = []

    for component in qualified_components:
        regions_in_group = [get_region(n) for n in component]
        centers = SkyCoord([r.center for r in regions_in_group])
        mean_center = SkyCoord(ra=centers.ra.mean(), dec=centers.dec.mean())

        # Average radius if CircleSkyRegion, otherwise use first
        try:
            mean_radius = np.mean([r.radius.to(u.arcsec).value for r in regions_in_group]) * u.arcsec
        except AttributeError:
            mean_radius = 1.0 * u.arcsec  # fallback default

        combined_region = CircleSkyRegion(center=mean_center, radius=mean_radius)
        combined_regions.append(combined_region)
        


    return combined_regions


def make_insignificant_catalogs(regions_b3, regions_b6, threshold=1.1e-5):
    """
    Match regions between B3 and B6 based on a threshold distance.
    
    Parameters:
    - regions_b3: Regions object for B3
    - regions_b6: Regions object for B6
    - threshold: Distance threshold for matching (default is 1.1e-5 degrees)
    
    Returns:
    - matched_regions: List of matched regions
    """
    
    matched_xsky = []
    matched_ysky = []
    matched_b6_idx = []
    matched_b3_idx = []
    b3only_xsky = []
    b3only_ysky = []
    for i, region_b3 in enumerate(regions_b3):
        center_b3 = region_b3.center
        
        for j, region_b6 in enumerate(regions_b6):
            center_b6 = region_b6.center
            
            if center_b3.separation(center_b6).value < threshold:
                matched_xsky.append(center_b6.ra.deg)
                matched_ysky.append(center_b6.dec.deg)
                matched_b6_idx.append(j)
                matched_b3_idx.append(i)
                break
    
        b3only_xsky.append(center_b3.ra.deg)
        b3only_ysky.append(center_b3.dec.deg)
    
    unmatched_b6_idx = [j for j in range(len(regions_b6)) if j not in matched_b6_idx]
    b6only_xsky = []
    b6only_ysky = []

    for j in unmatched_b6_idx:
        center_b6 = regions_b6[j].center
        b6only_xsky.append(center_b6.ra.deg)
        b6only_ysky.append(center_b6.dec.deg)       

    b3only_xsky = np.array(b3only_xsky)
    b3only_ysky = np.array(b3only_ysky)
    matched_xsky = np.array(matched_xsky)
    matched_ysky = np.array(matched_ysky)
    b6only_xsky = np.array(b6only_xsky)
    b6only_ysky = np.array(b6only_ysky)

    matched_tab = Table([matched_xsky, matched_ysky, np.ones(len(matched_xsky), dtype=bool), np.ones(len(matched_xsky), dtype=bool)], names=('ra', 'dec', 'isdetected_b6', 'isdetected_b3'))
    b3only_tab = Table([b3only_xsky, b3only_ysky, np.zeros(len(b3only_xsky), dtype=bool), np.ones(len(b3only_xsky), dtype=bool)], names=('ra', 'dec', 'isdetected_b6', 'isdetected_b3'))
    b6only_tab = Table([b6only_xsky, b6only_ysky, np.ones(len(b6only_xsky), dtype=bool), np.zeros(len(b6only_xsky), dtype=bool)], names=('ra', 'dec', 'isdetected_b6', 'isdetected_b3'))

    matched_regions = vstack([matched_tab, b3only_tab, b6only_tab])

    
    return matched_regions

def region_to_ds9_line(region):
    """Convert a SkyRegion to a DS9 format string line."""
    if region.__class__.__name__ == "CircleSkyRegion":
        ra = region.center.ra.deg
        dec = region.center.dec.deg
        radius = region.radius.to("deg").value
        return f"circle({ra:.8f},{dec:.8f},{radius:.8f})"
    
    elif region.__class__.__name__ == "EllipseSkyRegion":
        ra = region.center.ra.deg
        dec = region.center.dec.deg
        major = region.width.to("deg").value
        minor = region.height.to("deg").value
        angle = region.angle.to("deg").value
        return f"ellipse({ra:.8f},{dec:.8f},{major:.8f},{minor:.8f},{angle:.2f})"

    elif region.__class__.__name__ == "RectangleSkyRegion":
        ra = region.center.ra.deg
        dec = region.center.dec.deg
        width = region.width.to("deg").value
        height = region.height.to("deg").value
        angle = region.angle.to("deg").value
        return f"box({ra:.8f},{dec:.8f},{width:.8f},{height:.8f},{angle:.2f})"
    
    else:
        raise NotImplementedError(f"Region type {type(region)} not supported.")

def save_regions_to_ds9_manual(region_list, filename):
    with open(filename, "w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n")
        f.write("fk5\n")
        for region in region_list:
            try:
                line = region_to_ds9_line(region)
                f.write(line + "\n")
            except NotImplementedError as e:
                print(f"Skipping unsupported region: {e}")

if __name__ == "__main__":
    for reg in ['W51-E', 'W51-IRS2']:
        
        for band in ['B3', 'B6']:
            adam_file = f'tables/{reg}_{band}_adam_insignificant.reg'
            nazar_file = f'tables/{reg}_{band}_nazar_insignificant.reg'
            taehwa_file = f'tables/{reg}_{band}_taehwa_insignificant.reg'
            
            insign_regs = get_insignificant_regions(adam_file, nazar_file, taehwa_file)
            save_regions_to_ds9_manual(insign_regs, f'tables/{reg}_{band}_insignificant.reg')


        b3_insignificant = Regions.read(f'tables/{reg}_B3_insignificant.reg', format='ds9')
        b6_insignificant = Regions.read(f'tables/{reg}_B6_insignificant.reg', format='ds9')
        matched_insign = make_insignificant_catalogs(b3_insignificant, b6_insignificant)
        matched_insign.write(f'tables/{reg}_insignificant_catalog.fits', overwrite=True)

