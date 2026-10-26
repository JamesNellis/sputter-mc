import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from lattice import Lattice
from analysis import patch_fraction

# lattice
SIZE = 400
SEED = 0
SPACING = (1.0, 1.0)

# sputtering
EXTRA_ATOMS_CHOICES = (0, 1, 2)
THETA_MAX = 0.40
TOTAL_SPUTTERS = None
POINTS = 20

# analysis
R_VALUES = [2, 4, 6]

# output
OUT = "figs/run.png"
SHOW = True



here = Path(__file__).parent

atoms_per_ion = 1 + np.mean(EXTRA_ATOMS_CHOICES)
n_ions = TOTAL_SPUTTERS or round(THETA_MAX * SIZE**2 / atoms_per_ion)

lattice = Lattice(SIZE, seed=SEED, spacing=SPACING)
thetas = [0.0]
clean_fracs = {r: [1.0] for r in R_VALUES}

for batch in np.diff(np.linspace(0, n_ions, POINTS + 1).round().astype(int)):
    for _ in range(batch):
        lattice.sputter_one_ion(EXTRA_ATOMS_CHOICES)

    thetas.append(lattice.theta)

    for r, f in zip(R_VALUES, patch_fraction(lattice.damage, R_VALUES, lattice.spacing)):
        clean_fracs[r].append(f)

# print(f"{lattice.n_ions} ions, theta = {lattice.theta:.3f} ML, clean fraction = {lattice.clean_fraction:.3f}")


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
