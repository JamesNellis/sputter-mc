import matplotlib.pyplot as plt
from lattice import Lattice
from analysis import patch_fraction

lattice = Lattice(200, seed=0)
r_values = [2, 4, 6]

damaged_fractions = []
intensities = {r: [] for r in r_values}


total_sputters = 20000
points = 30

for i in range(points):
    for _ in range(total_sputters // points):
        lattice.sputter_one_ion()

    # x axis
    damaged_fractions.append(1 - lattice.clean_fraction)
    # eg [.8, .4, .2] - fractions for each r
    fracs = patch_fraction(lattice.damage, r_values)
    for r, f in zip(r_values, fracs):
        intensities[r].append(f)

for r in r_values:
    plt.plot(damaged_fractions, intensities[r], label=f"r={r}")

plt.xlabel("Damaged fraction")
plt.ylabel("Portion of cells which center a clean circle (radius = r)") # as in, portion of all cells which have a circle around them radius r, in which all cells are clean.
plt.title(f"{total_sputters} sputters")
plt.legend()
plt.savefig("figs/fig.png")
plt.show()
