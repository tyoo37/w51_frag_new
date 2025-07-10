import Paths.Paths as paths
from astropy.table import Table
Path = paths.filepaths()

w51e_matched_dendro = Table.read(Path.w51e_matched_dendro_catalog_new)
w51n_matched_dendro = Table.read(Path.w51n_matched_dendro_catalog_new)

