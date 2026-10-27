import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from lattice import Lattice
from analysis import patch_fraction

# lattice
SIZE = 400 # NxN lattice
SEED = 0
SPACING = (1.0, 1.0) # (row, col) spacing

# sputtering
EXTRA_ATOMS_CHOICES = (0, 1, 2)
THETA_MAX = 0.40 # Young's goes to .4 - sim aims for this much 
TOTAL_SPUTTERS = None # rn this is calculated at run based on theta_max
POINTS = 20 # how many samples on a curve - theta_max used to find total ions, then split based on points

# analysis
R_VALUES = [2, 4, 6] # tracking damage as defined by these r values

# output
OUT = "figs/run.png"
SHOW = True



here = Path(__file__).parent

# each ion removes 1 atom plus on average mean(EXTRA_ATOMS_CHOICES) neighbors
atoms_per_ion = 1 + np.mean(EXTRA_ATOMS_CHOICES)
n_ions = TOTAL_SPUTTERS or round(THETA_MAX * SIZE**2 / atoms_per_ion) # find total atoms, then divide by atoms/ion

lattice = Lattice(SIZE, seed=SEED, spacing=SPACING)
thetas = [0.0] # x values
clean_fracs = {r: [1.0] for r in R_VALUES} # coverage for each r value

# batch boundaries so the batches add up to exactly n_ions
for batch in np.diff(np.linspace(0, n_ions, POINTS + 1).round().astype(int)): # this just splits up n_ions in POINTS batches
    for _ in range(batch):
        lattice.sputter_one_ion(EXTRA_ATOMS_CHOICES)

    thetas.append(lattice.theta) # measure theta - x axis

    # measure clean fractions based on different rs - y axis
    for r, f in zip(R_VALUES, patch_fraction(lattice.damage, R_VALUES, lattice.spacing)):
        clean_fracs[r].append(f)

# print(f"{lattice.n_ions} ions, theta = {lattice.theta:.3f} ML, clean fraction = {lattice.clean_fraction:.3f}")


# draw lines w/ matplot
for r in R_VALUES:
    plt.plot(thetas, clean_fracs[r], label=f"r={r}")

plt.xlabel("θ (ML removed)")
plt.ylabel("fraction of sites centering a clean circle of radius r")
plt.title(f"{SIZE}x{SIZE}, spacing={tuple(round(float(s), 3) for s in SPACING)}, "
          f"extras={EXTRA_ATOMS_CHOICES}, {lattice.n_ions} ions")
plt.legend()

if OUT:
    (here / OUT).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(here / OUT)
    print(f"saved {here / OUT}")
if SHOW:
    plt.show()
