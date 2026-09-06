import numpy as np
from scipy.ndimage import distance_transform_edt


def patch_radii(damage):
    # for each clean cell finds the largest int r such that all cells within are clean
    # returns -1 for damaged cells

    size = damage.shape[0] # assumes square...
    clean = (damage == 0) # T/F grid same size as lattice

    # if no sputtering, just return a grid reporting rs of side of lattice
    if clean.all():
        return np.full_like(damage, size) 

    # no.tile takes the input - clean - and places it in the middle of 
    # a, in this case, 3x3 tiling of itself, so this solves the wrapping 
    # needed for the periodic boundaries used to imitate infinite x and y
    clean_tiled = np.tile(clean, (3, 3))

    # distance from each clean cell to nearest damaged cell. damaged cell: 0
    dist_tiled = distance_transform_edt(clean_tiled)

    # from the tiled copying, take out the center one
    dist = dist_tiled[size:2*size, size:2*size]

    # .where makes all damaged cells -1
    # np.ceil(dist).astype(int) - 1 finds the greatest int radius for which
    # all cells within the circle (inclusive) are clean.
    radii = np.where(clean, np.ceil(dist).astype(int) - 1, -1)
    return radii


def patch_fraction(damage, r_values):
    # fraction of cells with clean patch of at least radius r, for each r in r_values
    radii = patch_radii(damage)
    return [float(np.mean(radii >= r)) for r in r_values]
