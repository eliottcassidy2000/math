"""Referee: odd thetas theta_p^k + (p-1)K1 beyond the census.
(1) Z_p-mark obstruction (THM-4835 6(d)) for p in {3,5,7}, k in {2,3}.
(2) COUNT tier for theta_5^2 + 4K1 (A = S5 x S4): is |{T >= D}/A| equal to |Tour(X)/A| for some
    A-set X of size 10?  Necessary local test over odd-order classes (cycle types on X constrained
    by power maps); if no local profile matches, the count tier fails.  Run: python thetas.py"""
import itertools
from math import comb
from fractions import Fraction

from core import ordered_pair_orbits


def theta(p, k):
    u, w = 0, 1
    arcs, nxt = [], 2
    paths = []
    for _ in range(p):
        seq = [u] + list(range(nxt, nxt + k - 1)) + [w]
        paths.append(list(range(nxt, nxt + k - 1)))
        nxt += k - 1
        arcs += [(seq[i], seq[i + 1]) for i in range(len(seq) - 1)]
    isol = list(range(nxt, nxt + p - 1))
    return nxt + p - 1, arcs, paths, isol


def perm_from(n, path_perm, iso_perm, paths, isol):
    img = list(range(n))
    for i, P in enumerate(paths):
        for t, v in enumerate(P):
            img[v] = paths[path_perm[i]][t]
    for i, v in enumerate(isol):
        img[v] = isol[iso_perm[i]]
    return tuple(img)


def orbits_pairs_of_cycle_type(ct, npts):
    """number of <g>-orbits on unordered pairs for g with cycle type ct (dict len->count)."""
    img = []
    base = 0
    for L, c in sorted(ct.items()):
        for _ in range(c):
            img += [base + (t + 1) % L for t in range(L)]
            base += L
    assert base == npts
    o, rev = ordered_pair_orbits([tuple(img)], npts, frozenset())
    return o  # None if reversing (even cycles) -- not used for odd


def zp_mark(p, k):
    n, arcs, paths, isol = theta(p, k)
    E = {tuple(sorted(a)) for a in arcs}
    rho = perm_from(n, [(i + 1) % p for i in range(p)], list(range(p - 1)), paths, isol)
    o_pi, _ = ordered_pair_orbits([rho], n, E)
    m = n - 1
    possible = []
    for x in range(0, m + 1):
        if (m - x) % p:
            continue
        o_x = (comb(m, 2) + (p - 1) * comb(x, 2)) // p
        possible.append((x, o_x))
    match = [x for x, o in possible if o == o_pi]
    fixedV = sum(1 for v in range(n) if rho[v] == v)
    print(f"theta_{p}^{k}+{p - 1}K1 (n={n}): Z_{p}-orbits on Pi = {o_pi} (|V^rho|={fixedV}); "
          f"X options (x, orbits) = {possible}; matching x: {match}")


def count_test_5_2():
    p, k = 5, 2
    n, arcs, paths, isol = theta(p, k)
    E = {tuple(sorted(a)) for a in arcs}
    m = n - 1
    # odd classes of S5 x S4: representatives (path perm, isolate perm) and class sizes
    c3_5 = (1, 2, 0, 3, 4)
    c5_5 = (1, 2, 3, 4, 0)
    id5 = (0, 1, 2, 3, 4)
    c3_4 = (1, 2, 0, 3)
    id4 = (0, 1, 2, 3)
    classes = {"11": (id5, id4, 1), "31": (c3_5, id4, 20), "51": (c5_5, id4, 24),
               "13": (id5, c3_4, 8), "33": (c3_5, c3_4, 160), "53": (c5_5, c3_4, 192)}
    order = sum(c[2] for c in classes.values())
    assert order == 405
    lhs = 0
    opi = {}
    for name, (pp, ip, sz) in classes.items():
        g = perm_from(n, pp, ip, paths, isol)
        o, rev = ordered_pair_orbits([g], n, E)
        assert not rev
        opi[name] = o
        lhs += sz * 2 ** o
    count = Fraction(lhs, 120 * 24)
    print(f"theta_5^2+4K1: odd-class orbit numbers on Pi {opi}; marked count = {count}")
    # an X on which some involution of A acts trivially has sum >= 2*2^45 > lhs, and if no
    # involution acts trivially every even-order element has an even cycle on X (fixes nothing);
    # so only odd-order classes of A contribute on the X side for a possible match.
    print(f"   lhs - 2^45 = {lhs - 2 ** 45} < 2^45: {lhs - 2 ** 45 < 2 ** 45}")
    matches = []
    for xa in (1, 4, 7, 10):          # class 31 (order 3)
        for xb in (1, 4, 7, 10):      # class 13
            for xc in (1, 4, 7, 10):  # class 33
                for xd in (0, 5, 10):  # class 51 (order 5)
                    c1 = xd + xb - m  # class 53 (order 15): fixed points
                    if c1 < 0 or (xd - c1) % 3 or (xb - c1) % 5:
                        continue
                    c3, c5 = (xd - c1) // 3, (xb - c1) // 5
                    if c1 + 3 * c3 + 5 * c5 != m:
                        continue
                    ox = {"11": comb(m, 2),
                          "31": orbits_pairs_of_cycle_type({1: xa, 3: (m - xa) // 3}, m),
                          "13": orbits_pairs_of_cycle_type({1: xb, 3: (m - xb) // 3}, m),
                          "33": orbits_pairs_of_cycle_type({1: xc, 3: (m - xc) // 3}, m),
                          "51": orbits_pairs_of_cycle_type({1: xd, 5: (m - xd) // 5}, m),
                          "53": orbits_pairs_of_cycle_type({k_: v for k_, v in {1: c1, 3: c3, 5: c5}.items() if v}, m)}
                    rhs = sum(classes[c][2] * 2 ** ox[c] for c in classes)
                    if rhs == lhs:
                        matches.append((xa, xb, xc, xd, c1))
    print(f"   local odd-class profiles of 10-point actions giving the same count: {matches}")


if __name__ == "__main__":
    for p in (3, 5, 7):
        for k in (2, 3):
            zp_mark(p, k)
    count_test_5_2()
