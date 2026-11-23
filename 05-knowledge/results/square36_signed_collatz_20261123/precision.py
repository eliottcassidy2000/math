"""Exact binary stopping precisions for one flip from height 199."""
from math import isqrt

from physical import path_count, strong_catalogue


def main():
    h = 199
    cores = strong_catalogue()[h]
    # The only move entirely in exterior singletons with gain 2 mod 4 has
    # one intervening singleton and multiplies the whole height by three.
    readings = {3*h}
    for core in cores:
        readings.update(core.flips)
        for marked, full in ((core.starts, core.end), (core.ends, core.start)):
            for mass, old in zip(marked, full):
                readings.add(h+mass-old)
                readings.add(h+2*mass)
    readings = sorted(x for x in readings if (x-h) % 4 == 2)
    print("Strong cores", len(cores), "opposite-type one-flip heights", len(readings))
    assert len(cores) == 39 and len(readings) == 347
    square_neighbors = [x for x in readings if isqrt(h+x)**2 == h+x]
    assert square_neighbors == [57, 125, 201, 285, 377, 585, 701, 825, 1097, 1245, 2505, 3897, 4157, 6857]
    print("All square-sum one-flip neighbors", square_neighbors)
    for target in (201, 701):
        assert target in readings
        k = 2
        while any((x-target) % (1 << k) == 0 for x in readings if x != target):
            k += 1
        last = [x for x in readings if x != target and (x-target) % (1 << (k-1)) == 0]
        assert k == 10 and last == ({201: [713, 3785, 6857], 701: [189, 2237]}[target])
        square_k = max(((x-target) & -(x-target)).bit_length()
                       for x in square_neighbors if x != target)
        assert square_k == (10 if target == 201 else 8)
        print("Target", target, "least binary precision", k, "last aliases", last)

    rows = [510, 508, 128, 228, 172, 388, 436, 256, 28]
    for first, other, middle in [((5,6),(0,3),57), ((0,3),(5,6),141)]:
        host = rows.copy()
        assert path_count(host) == 43
        for (u, v), height in [(first, middle), (other, 199)]:
            host[u] ^= 1 << v
            host[v] ^= 1 << u
            assert path_count(host) == height
    lower = (1 << 60)*((1 << 61)-1)-1
    assert lower+3 < 1 << 121 and 5**156 <= lower**3 < 5**157


if __name__ == "__main__":
    main()
