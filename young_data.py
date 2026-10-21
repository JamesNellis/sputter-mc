import matplotlib.pyplot as plt
import numpy as np

# first column is ML removed, rest are f for each r
lines = open("data/MC_Data.txt").read().splitlines()
headers = [h.strip("'") for h in lines[0].split("\t")]  # some headers have quotes,
data = np.array([[float(x) for x in line.split("\t")] for line in lines[1:] if line.strip()])

theta = data[:, 0]
for r in ["r2", "r3", "r4", "r5", "r6"]:
    plt.plot(theta, data[:, headers.index(r)], "o-", label=r)

plt.xlabel("ML removed (Young's theta)")
plt.ylabel("fraction of sites centering a clean circle of radius r")
plt.title("Young's MC data")
plt.legend()
plt.savefig("figs/young.png")
plt.show()
