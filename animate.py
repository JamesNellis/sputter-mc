import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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
FRAMES = 80                     # one frame per batch of ions

# analysis
R_VALUES = [2, 4, 6]

# animation
VIEW = 120                      # show a VIEW x VIEW corner of the lattice
FRAME_MS = 80                   # time per frame
HOLD_MS = 1500                  # pause on the last frame before looping

# output
OUT = "figs/sputter.gif"
SHOW = True


# run the whole simulation
atoms_per_ion = 1 + np.mean(EXTRA_ATOMS_CHOICES)
n_ions = round(THETA_MAX * SIZE**2 / atoms_per_ion)

lattice = Lattice(SIZE, seed=SEED, spacing=SPACING)
thetas = [0.0]
clean_fracs = {r: [1.0] for r in R_VALUES}
frames = [lattice.damage[:VIEW, :VIEW].copy()]  # damage in the corner, one per frame

print(f"simulating {n_ions:,} ions over {FRAMES} frames...")
for batch in np.diff(np.linspace(0, n_ions, FRAMES + 1).round().astype(int)):
    for _ in range(batch):
        lattice.sputter_one_ion(EXTRA_ATOMS_CHOICES)
    thetas.append(lattice.theta)
    for r, f in zip(R_VALUES, patch_fraction(lattice.damage, R_VALUES, lattice.spacing)):
        clean_fracs[r].append(f)
    frames.append(lattice.damage[:VIEW, :VIEW].copy())


fig, (ax_lat, ax_frac) = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")

# left: the lattice
sy, sx = SPACING
top = frames[-1].max()
# damage is a whole number of atoms, so one color per count
im = ax_lat.imshow(frames[0], cmap=plt.get_cmap("viridis", top + 1), vmin=-0.5, vmax=top + 0.5,
                   extent=[0, VIEW * sx, VIEW * sy, 0], interpolation="nearest")
fig.colorbar(im, ax=ax_lat, label="atoms removed", ticks=range(top + 1))
ax_lat.set_xlabel("x (a)")
ax_lat.set_ylabel("y (a)")
ax_lat.set_title(f"θ = {thetas[0]:.3f} ML")

# right: clean fraction vs theta, full curves with a line at the current frame's theta
for r in R_VALUES:
    ax_frac.plot(thetas, clean_fracs[r], label=f"r={r}")
now = ax_frac.axvline(thetas[0], color="k", lw=1)
ax_frac.set_xlim(0, thetas[-1])
ax_frac.set_ylim(0, 1.02)
ax_frac.set_xlabel("θ (ML removed)")
ax_frac.set_ylabel("fraction of sites centering a clean circle of radius r")
ax_frac.legend(loc="upper right")  # fixed spot, the default "best" moves it out of the line's way every frame


def draw(i):
    im.set_data(frames[i])
    ax_lat.set_title(f"θ = {thetas[i]:.3f} ML")
    now.set_xdata([thetas[i], thetas[i]])
    return [im, now]


here = Path(__file__).parent

if OUT:
    # a gif can't pause before looping, so repeat the last frame for HOLD_MS instead
    order = list(range(len(frames))) + [len(frames) - 1] * (HOLD_MS // FRAME_MS)
    FuncAnimation(fig, draw, frames=order).save(here / OUT, writer="pillow", fps=1000 / FRAME_MS, dpi=100)
    print(f"saved {here / OUT}")
if SHOW:
    anim = FuncAnimation(fig, draw, frames=len(frames), interval=FRAME_MS, repeat_delay=HOLD_MS)
    plt.show()
