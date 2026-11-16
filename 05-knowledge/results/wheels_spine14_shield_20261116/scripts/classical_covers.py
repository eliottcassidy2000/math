#!/usr/bin/env python3
"""classical_covers.py: (1) Erdos's 1950 covering system {0 (2), 0 (3), 1 (4), 3 (8), 7 (12), 23 (24)} covers Z/24 and is the
first pattern of CLOCKS.md for the Riesel wheel (orders 2,4,3,12,8,24 with phases 0,1,0,7,3,23); (2) the Brier parity bound
sum_keepers 1/d + 2 sum_flippers(d>2) 1/d + (1/2) sum_odd 1/d >= 1 (with 3 in the set) evaluated on W(24), W(36), W(64) and on
the Filaseta-Finch-Kozek Brier covering set of period 720."""
from sympy import n_order
erdos = {2: 0, 3: 0, 4: 1, 8: 3, 12: 7, 24: 23}
covered = all(any(x % d == r for d, r in erdos.items()) for x in range(24))
print("Erdos 1950 system covers Z/24:", covered, "; primes of the moduli:", {d: p for d, p in zip((2, 3, 4, 8, 12, 24), (3, 7, 5, 17, 13, 241))})
clocks_order = [2, 4, 3, 12, 8, 24]; clocks_phase = [0, 1, 0, 7, 3, 23]
print("CLOCKS.md pattern (0,1,0,7,3,23) in order (2,4,3,12,8,24) equals Erdos's system:", {d: r for d, r in zip(clocks_order, clocks_phase)} == erdos)
def brier_bound(primes):
    d = {p: n_order(2, p) for p in primes}
    keep = sum(1 / d[p] for p in primes if d[p] % 4 == 0)
    flip = sum(1 / d[p] for p in primes if d[p] % 4 == 2 and d[p] > 2)
    odd = sum(1 / d[p] for p in primes if d[p] % 2 == 1)
    return keep + 2 * flip + odd / 2, (3 in primes)
for name, S in (("W(24)", [3, 5, 7, 13, 17, 241]), ("W(36)", [3, 5, 7, 13, 19, 37, 73, 109]), ("W(64)", [3, 5, 17, 257, 641, 65537, 6700417]),
                ("FFK Brier set (period 720)", [3, 5, 7, 11, 13, 17, 19, 31, 37, 41, 61, 73, 97, 109, 151, 241, 257, 331, 673, 1321])):
    val, has3 = brier_bound(S)
    print(f"{name}: Brier parity bound = {val:.3f} (needs >= 1; 3 in set: {has3}) -> {'possible' if val >= 1 else 'no Brier number can have a covering set inside it'}")
