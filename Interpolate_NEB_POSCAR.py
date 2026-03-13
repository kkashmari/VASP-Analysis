import os
import argparse
import numpy as np


def read_poscar(path):
    with open(path, "r") as f:
        lines = [l.strip() for l in f.readlines()]

    comment = lines[0]
    scale = float(lines[1])
    lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))
    coord_type = lines[7].strip()

    natoms = sum(counts)
    positions = np.array(
        [[float(x) for x in lines[8 + i].split()] for i in range(natoms)]
    )

    return comment, scale, lattice, elements, counts, coord_type, positions


def write_poscar(path, comment, scale, lattice, elements, counts, coord_type, positions):
    with open(path, "w") as f:
        f.write(comment + "\n")
        f.write(f"{scale:16.10f}\n")
        for v in lattice:
            f.write(" ".join(f"{x:16.10f}" for x in v) + "\n")
        f.write(" ".join(elements) + "\n")
        f.write(" ".join(map(str, counts)) + "\n")
        f.write(coord_type + "\n")
        for p in positions:
            f.write(" ".join(f"{x:16.10f}" for x in p) + "\n")


def pbc_delta(delta):
    """Apply minimum image convention in fractional coordinates"""
    delta[delta > 0.5] -= 1.0
    delta[delta < -0.5] += 1.0
    return delta


def interpolate(start, end, t):
    return start + t * pbc_delta(end - start)


def main():
    parser = argparse.ArgumentParser(description="Generate NEB interpolated POSCARs")
    parser.add_argument("-n", "--nimages", type=int, required=True,
                        help="Total number of NEB images (including endpoints)")
    args = parser.parse_args()

    n = args.nimages

    c0, s0, lat0, elems, nums, ctype, pos0 = read_poscar("00/POSCAR")
    c1, s1, lat1, _, _, _, pos1 = read_poscar(f"{n-1:02d}/POSCAR")

    for i in range(n):
        t = i / (n - 1)

        os.makedirs(f"{i:02d}", exist_ok=True)

        scale = (1 - t) * s0 + t * s1
        lattice = (1 - t) * lat0 + t * lat1
        positions = interpolate(pos0, pos1, t)

        write_poscar(
            f"{i:02d}/POSCAR",
            c0,
            scale,
            lattice,
            elems,
            nums,
            ctype,
            positions,
        )


if __name__ == "__main__":
    main()
