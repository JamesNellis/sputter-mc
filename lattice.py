import numpy as np

class Lattice:
    def __init__(self, size, seed):
        self.size = size
        self.seed = seed
        self.damage = np.zeros((size, size), dtype=int)
        self.rng = np.random.default_rng(seed)
        self.n_ions = 0

    def random_cell(self):
        row = self.rng.integers(0, self.size)
        col = self.rng.integers(0, self.size)
        return (row, col)

    def remove_atom(self, cell):
        self.damage[cell] += 1

    def neighbors(self, cell):
        row, col = cell

        # this version wraps: 0,0 will return 4 neighbors still
        return [
            ((row - 1) % self.size, col),
            ((row + 1) % self.size, col),
            (row, (col - 1) % self.size),
            (row, (col + 1) % self.size),
        ]

    def sputter_one_ion(self, extra_atoms_choices=(0, 1, 2)):
    # removes an atom from the one random cell, then also removes an atom 
        # from n_extra_atoms neighbors, n_extra_atoms being a random choice of extra_atoms_choices
        cell = self.random_cell()
        self.remove_atom(cell)
        n_extra_atoms = self.rng.choice(extra_atoms_choices)
        if n_extra_atoms > 0:
            neighbor_list = self.neighbors(cell)
            indices = self.rng.choice(4, size=n_extra_atoms, replace=False)
            for idx in indices:
                self.remove_atom(neighbor_list[idx])
        self.n_ions += 1


    @property
    def clean_fraction(self):
        # fraction of surface "clean" or 0
        return (self.damage == 0).mean()


if __name__ == "__main__":
    lattice = Lattice(200, seed=0)
    for _ in range(10_000):
        lattice.sputter_one_ion()
    print("damage.sum():", lattice.damage.sum())
    print("(damage > 0).sum():", (lattice.damage > 0).sum())
    print("damage.max():", lattice.damage.max())
    print("n_ions:", lattice.n_ions)
    print("clean portion:", lattice.clean_fraction)
