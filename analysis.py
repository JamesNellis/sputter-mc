import numpy as np
from scipy.ndimage import distance_transform_edt


def patch_radii(damage, spacing=(1.0, 1.0)):
    # for each clean cell finds the distance (in units of a) to the nearest damaged cell,
    # so every cell strictly closer than that is clean
    # returns -1 for damaged cells

    size = damage.shape[0] # assumes square...
    clean = (damage == 0) # T/F grid same size as lattice

    # if no sputtering, return array of inf
    if clean.all():
        return np.full(damage.shape, np.inf)

    # no.tile takes the input - clean - and places it in the middle of 
    # a, in this case, 3x3 tiling of itself, so this solves the wrapping 
    # needed for the periodic boundaries used to imitate infinite x and y

    # could have also jyst made the lattice bigger than needed, and cut out the frame before any analysis.
    clean_tiled = np.tile(clean, (3, 3))

    # distance from each clean cell to nearest damaged cell. damaged cell: 0
    # sampling= changes the coefficients that are used in pythagorean theorem

    dist_tiled = distance_transform_edt(clean_tiled, sampling=spacing)

    # from the tiled copying, take out the center one
    dist = dist_tiled[size:2*size, size:2*size]

    # .where makes all damaged cells -1
    radii = np.where(clean, dist, -1.0)
    return radii


def patch_fraction(damage, r_values, spacing=(1.0, 1.0)):
    # fraction of cells with clean patch of at least radius r, for each r in r_values
    radii = patch_radii(damage, spacing)
    return [float(np.mean(radii > r)) for r in r_values]
