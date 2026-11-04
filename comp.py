import numpy as np
import matplotlib.pyplot as plt
from lattice import Lattice
from analysis import patch_fraction

# adjust for problem
adj_scale = 5 # 1 = no adjustment, just using data as given

# settings
SIZE = 500
SEED = 0
SPACING = (1.0, 1/np.sqrt(2))  # from Ni(110) (?)
EXTRA_ATOMS_CHOICES = (0, 1, 2)
R_VALUES = [2, 3, 4, 5, 6] # analysis
OUT = "figs/compare.png" if adj_scale == 1 else f"figs/compare_adj{adj_scale}.png"

# Young's data: first row is r values, first column is his theta (damage.sum() / size^2)
d = np.genfromtxt("data/young.csv", delimiter=",")
young_r, young_theta, young_f = d[0, 1:], d[1:, 0], d[1:, 1:]

# simulation, sputtered until max theta
lattice = Lattice(SIZE, seed=SEED, spacing=SPACING)
thetas = [0.0]
mine = {r: [1.0] for r in R_VALUES}
for target in np.linspace(0, young_theta.max() / adj_scale, 40)[1:]: # pick 40 evenly spaced points between 0 and youngs max theta
    while lattice.theta < target:
        lattice.sputter_one_ion(EXTRA_ATOMS_CHOICES)
    thetas.append(lattice.theta)
    for r, f in zip(R_VALUES, patch_fraction(lattice.damage, R_VALUES, lattice.spacing)):
        mine[r].append(f)

# plot: lines = sim, circles = young data
for i, r in enumerate(R_VALUES):
    plt.plot(thetas, mine[r], color=f"C{i}", label=f"r={r}")
    plt.plot(young_theta / adj_scale, young_f[:, np.argmin(abs(young_r - r))], "o", mfc="none", color=f"C{i}")

plt.plot([], [], "k-", label="lines: simulation")
plt.plot([], [], "ko", mfc="none", label="circles: Young data")
plt.xlabel("θ (ML removed)")
plt.ylabel("fraction of sites centering a clean circle of radius r")
plt.title(f"Simulation vs Young's data, Ni(110), {SIZE}x{SIZE}")
plt.legend(fontsize=8)
plt.savefig(OUT)
plt.show()