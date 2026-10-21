import numpy as np
import matplotlib.pyplot as plt
from lattice import Lattice
from analysis import patch_fraction

lattice = Lattice(400, seed=0, spacing=(1.0, 1/np.sqrt(2)))  # Ni(110)
r_values = [2, 4, 6]

thetas = []
intensities = {r: [] for r in r_values}


total_sputters = 33000  # ~0.41 ML on 400x400 at ~2 atoms per ion, matching Young's range
points = 20

for i in range(points):
    for _ in range(total_sputters // points):
        lattice.sputter_one_ion()

    # x axis
    thetas.append(lattice.theta)
    # eg [.8, .4, .2] - fractions for each r
    fracs = patch_fraction(lattice.damage, r_values, lattice.spacing)
    for r, f in zip(r_values, fracs):
        intensities[r].append(f)

for r in r_values:
    plt.plot(thetas, intensities[r], label=f"r={r}")

plt.xlabel("damage.sum() / number of sites")
plt.ylabel("Portion of cells which center a clean circle (radius = r)") # as in, portion of all cells which have a circle around them radius r, in which all cells are clean.
plt.title(f"{total_sputters} sputters")
plt.legend()
plt.savefig("figs/fig.png")
plt.show()
