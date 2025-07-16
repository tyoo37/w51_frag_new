from regions import Regions
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

w51e_b6_adam_unselected = Regions.read('tables/w51e_b6_adam_insignificant.reg', format='ds9')
nazar_unselected = Regions.read('tables/w51e_b6_nazar_insignificant.reg', format='ds9')
taehwa_unselected = Regions.read('tables/w51e_b6_taehwa_insignificant.reg', format='ds9')


def get_insignificant_regions(region_a_file, region_b_file, region_c_file, match_radius=1.1e-5):
    regions_A = Regions.read(region_a_file)
    regions_B = Regions.read(region_b_file)
    regions_C = Regions.read(region_c_file)

    all_regions = {'A': regions_A, 'B': regions_B, 'C': regions_C}

    # Define a match radius for overlap (e.g., 1 arcsec)

    # Function to get center coordinate
    def get_center(region):
        return region.center

    # Store matches: key = region, value = set of other files it overlaps with
    matches = defaultdict(set)

    # Compare across files
    file_names = list(all_regions.keys())

    for i, name1 in enumerate(file_names):
        for j, name2 in enumerate(file_names):
            if i >= j:
                continue  # Avoid self-comparison and double-comparing
            regs1 = all_regions[name1]
            regs2 = all_regions[name2]

            for r1 in regs1:
                c1 = get_center(r1)
                for r2 in regs2:
                    c2 = get_center(r2)
                    sep = c1.separation(c2)
                    if sep < match_radius:
                        matches[r1].add(name2)
                        matches[r2].add(name1)

    # Keep only regions that match with >1 other files
    final_regions = [r for r, matched_files in matches.items() if len(matched_files) >= 2]


    return final_regions


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
            
            if center_b3.separation(center_b6) < threshold:
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


if __name__ == "__main__":
    for reg in ['w51e', 'w51n']:
        
        for band in ['b3', 'b6']:
            adam_file = f'tables/{reg}_{band}_adam_insignificant.reg'
            nazar_file = f'tables/{reg}_{band}_nazar_insignificant.reg'
            taehwa_file = f'tables/{reg}_{band}_taehwa_insignificant.reg'
            
            insign_regs = get_insignificant_regions(adam_file, nazar_file, taehwa_file)
            insign_regs.write(f'tables/{reg}_{band}_insignificant.reg', format='ds9', overwrite=True)


        b3_insignificant = Regions.read(f'tables/{reg}_b3_insignificant.reg', format='ds9')
        b6_insignificant = Regions.read(f'tables/{reg}_b6_insignificant.reg', format='ds9')
        matched_insign = make_insignificant_catalogs(b3_insignificant, b6_insignificant)
        matched_insign.write(f'tables/{reg}_insignificant_catalog.fits', overwrite=True)

