#!/usr/bin/env python3
"""S91 (THM-5129): the complete supply of one-vertex sources of small offset, at every order.

A one-vertex source is a rooted tournament (B, y) with offset D = H(B) and slope U = U_B(y) = number of ordered covers
(P, R) of B by vertex-disjoint paths with R starting at y (P may be empty).  The closed tournament B => TT_t => v with
v -> y has H = D + 2^t U (t >= 1; near-join of HYP-9453, tail law of THM-5126 (A)).

Product law (KNOWN: Session216 star theorem, THM-5130 (F2); re-derived in THM-5129).  If B = C_1 => ... => C_m is the
decomposition into strong components and y lies in C_k, then
        D = prod_i H(C_i),      U = prod_{i<k} H(C_i) * U_{C_k}(y) * prod_{i>k} Pi_2(C_i),
Pi_2(C) = number of ordered pairs (P, R) of disjoint paths covering C, either possibly empty (Pi_2(point) = 2).
Hence the slopes of offset D are  Slopes_0(D) * 2^j (j >= 0),  with
        Slopes_0(D) = { d1 * u * w :  D = d1 * h * d3,  d1 a product of strong path counts,
                        u in Root(h) = {U_C(y) : C strong, H(C) = h} (Root(1) = {1}),
                        w in Suf(d3) = {prod Pi_2(C_i) : d3 = prod H(C_i), C_i strong of order >= 3} }.
A strong tournament of order n has at least f(n) Hamiltonian paths, f(3..9) = 3, 5, 9, 15, 25, 45, 75 (exhaustive
here for n <= 8; f(9) = 75 is checked on the list codes9.npy when present, and is Busch's theorem (2006) in general),
and f is nondecreasing.  So every strong component of a base with D <= 73 has order <= 8, and the lists below are
complete at every order of B.

Output (FINITE-EXACT + the product law):
  (1) the strong data of orders 3..8: numbers of strong classes, least H, numbers of distinct triples;
  (2) cross-check: for every D <= 73 the slopes of the complete rooted census of order <= 9 (rooted_le9.npz, when
      present) are contained in the product-law set, with equality when every reduced base has order <= 9;
  (3) the odd parts of the slopes at the Fermat-type offsets 3, 5, 9, 17, 33, 65 with their least 2-exponent j:
      the d-grid row of level s is supplied for all r >= 1 by a one-vertex source of offset 2^s + 1 iff j <= s;
  (4) the two-step / near-join values D + 2^t U (t >= 1) below 225: exactly the odd numbers other than
      1, 7, 11, 13, 21, 47, 67 (so these seven are missing at every order of the base);
  (5) the window sets of HYP-9453 (2): for every attained h <= 73 the number of slopes in [h, 3h/2].
Exit status 1 on any failed check.
"""
import os
import sys
from collections import defaultdict

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
DMAX = 73


@njit(cache=True)
def strong_data(n, reps, res):
    """res[r] = (H, Pi_2, strong flag, U(y) for y = 0..n-1) for every class representative of order n."""
    full = (1 << n) - 1
    g = np.zeros((1 << n, n), np.int64)
    f = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    out = np.zeros(n, np.int64)
    for r in range(reps.shape[0]):
        code = reps[r]
        for i in range(n):
            out[i] = 0
        for j in range(1, n):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    out[i] |= (1 << j)
                else:
                    out[j] |= (1 << i)
        hS[0] = 1
        for S in range(1, full + 1):
            tot = 0
            for v in range(n):
                f[S, v] = 0
                g[S, v] = 0
                if ((S >> v) & 1) == 0:
                    continue
                R = S ^ (1 << v)
                if R == 0:
                    f[S, v] = 1
                    g[S, v] = 1
                    tot += 1
                    continue
                a = 0
                b = 0
                for w in range(n):
                    if (R >> w) & 1:
                        if (out[w] >> v) & 1:
                            a += f[R, w]
                        else:
                            b += g[R, w]
                f[S, v] = a
                g[S, v] = b
                tot += a
            hS[S] = tot
        res[r, 0] = hS[full]
        pi2 = 0
        for S in range(0, full + 1):
            pi2 += hS[S] * hS[full ^ S]
        res[r, 1] = pi2
        # strong iff vertex 0 reaches everything and everything reaches vertex 0
        seen = 1
        frontier = 1
        while frontier:
            nf = 0
            for v in range(n):
                if (frontier >> v) & 1:
                    nf |= out[v]
            frontier = nf & ~seen
            seen |= nf
        fwd = seen == full
        seen = 1
        frontier = 1
        while frontier:
            nf = 0
            for v in range(n):
                if ((seen >> v) & 1) == 0:
                    if out[v] & frontier:
                        nf |= (1 << v)
            frontier = nf
            seen |= nf
        res[r, 2] = 1 if (fwd and seen == full) else 0
        for y in range(n):
            u = 0
            for S in range(1, full + 1):
                if (S >> y) & 1:
                    u += g[S, y] * hS[full ^ S]
            res[r, 3 + y] = u


def odd_part(u):
    j = (u & -u).bit_length() - 1
    return u >> j, j


def main():
    ok = True
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    root = defaultdict(set)      # h -> {U_C(y)}
    ncls = defaultdict(int)      # h -> number of strong classes with H = h
    pi2 = defaultdict(set)       # h -> {Pi_2(C)}
    least = {}
    print("(1) strong tournaments: order, classes, strong classes, least H, distinct (H, Pi_2), distinct (H, U)")
    for n in range(3, 9):
        reps = z["n%d" % n].astype(np.int64)
        res = np.zeros((len(reps), 3 + n), np.int64)
        strong_data(n, reps, res)
        st = res[res[:, 2] == 1]
        least[n] = int(st[:, 0].min())
        pairs_pi = set()
        pairs_u = set()
        for row in st:
            h = int(row[0])
            ncls[h] += 1
            pi2[h].add(int(row[1]))
            pairs_pi.add((h, int(row[1])))
            for y in range(n):
                root[h].add(int(row[3 + y]))
                pairs_u.add((h, int(row[3 + y])))
        print(f"    n = {n}: {len(reps):5d} classes, {len(st):5d} strong, least H = {least[n]:3d}, "
              f"{len(pairs_pi):5d} pairs (H, Pi_2), {len(pairs_u):6d} pairs (H, U)")
    ok &= [least[n] for n in range(3, 9)] == [3, 5, 9, 15, 25, 45]
    c9 = os.path.join(HERE, "codes9.npy")
    if os.path.exists(c9):
        codes9 = np.load(c9).astype(np.int64)
        res = np.zeros((len(codes9), 3 + 9), np.int64)
        strong_data(9, codes9, res)
        f9 = int(res[res[:, 2] == 1][:, 0].min())
        print(f"    n = 9: least H of a strong tournament = {f9} (list codes9.npy, every class of order 9)")
        ok &= f9 == 75
    else:
        print("    n = 9: codes9.npy not present; f(9) = 75 is taken from Busch (2006) / census.py")
    strong_vals = sorted(h for h in root if h <= DMAX)
    # Pre: products of strong values (the attained numbers <= DMAX); Suf: products of Pi_2 over ordered factorisations
    pre = {1}
    suf = {1: {1}}
    for d in range(3, DMAX + 1, 2):
        s = set()
        for h in strong_vals:
            if d % h == 0 and (d // h) in suf:
                for w in suf[d // h]:
                    for p in pi2[h]:
                        s.add(p * w)
        if s:
            suf[d] = s
            pre.add(d)
    print(f"    attained odd numbers <= {DMAX}: all odd numbers except {sorted(set(range(1, DMAX + 1, 2)) - pre)}")
    ok &= sorted(set(range(1, DMAX + 1, 2)) - pre) == [7, 21]
    slopes0 = {}
    for D in sorted(pre):
        s = set()
        for d1 in pre:
            if D % d1:
                continue
            rest = D // d1
            for h in [1] + strong_vals:
                if rest % h:
                    continue
                d3 = rest // h
                if d3 not in suf:
                    continue
                us = {1} if h == 1 else root[h]
                for u in us:
                    for w in suf[d3]:
                        s.add(d1 * u * w)
        slopes0[D] = s
    prim = {}
    for D, s in slopes0.items():
        m = {}
        for u in s:
            dl, j = odd_part(u)
            if dl not in m or j < m[dl]:
                m[dl] = j
        prim[D] = m

    def is_slope(D, u):
        if D not in prim:
            return False
        dl, j = odd_part(u)
        return dl in prim[D] and j >= prim[D][dl]

    # (2) cross-check with the complete rooted census of order <= 9
    cz = os.path.join(HERE, "rooted_le9.npz")
    if os.path.exists(cz):
        zz = np.load(cz)
        rows = np.concatenate([zz[f"m{m}"] for m in range(1, 10)], axis=0)
        cen = defaultdict(set)
        sel = rows[:, 0] <= DMAX
        for d, u in zip(rows[sel, 0].tolist(), rows[sel, 4].tolist()):
            cen[d].add(u)
        bad = [(d, u) for d in cen for u in cen[d] if not is_slope(d, u)]
        missing = {}
        for D in sorted(prim):
            cm = {}
            for u in cen.get(D, ()):
                dl, j = odd_part(u)
                if dl not in cm or j < cm[dl]:
                    cm[dl] = j
            diff = sorted(dl for dl in prim[D] if dl not in cm or cm[dl] > prim[D][dl])
            if diff:
                missing[D] = len(diff)
                print(f"    offset {D}: odd parts (with least exponent) not reached by a base of order <= 9: "
                      f"{[(dl, prim[D][dl]) for dl in diff]}")
        print(f"(2) census of order <= 9: {sum(len(v) for v in cen.values())} pairs with D <= {DMAX}; "
              f"not in the product-law set: {len(bad)}; offsets where the product law has odd parts (or smaller "
              f"exponents) beyond order 9: {missing}")
        ok &= len(bad) == 0
        ok &= set(cen) == set(prim)
    else:
        print("(2) rooted_le9.npz not present: cross-check skipped (regenerate with census.py)")
    # (3) Fermat-type offsets
    print("(3) odd parts d of the slopes at offset 2^s + 1, with the least 2-exponent j (d:j); '*' = j > s (finite head)")
    for s in range(1, 7):
        D = (1 << s) + 1
        m = prim[D]
        cells = [f"{dl}:{j}{'*' if j > s else ''}" for dl, j in sorted(m.items())]
        full_rows = sorted(dl for dl, j in m.items() if j <= s)
        print(f"    s = {s}, D = {D}: {len(m)} odd parts, {len(full_rows)} with j <= s")
        line = "        "
        for c in cells:
            if len(line) + len(c) > 118:
                print(line)
                line = "        "
            line += c + " "
        print(line)
    ok &= prim[3] == {3: 0}
    ok &= prim[5] == {1: 5, 3: 1, 5: 0, 7: 0, 9: 0}
    ok &= prim[17] == {5: 2, 7: 2, 9: 1, 11: 2, 17: 0, 21: 0, 31: 0, 33: 0, 55: 0, 61: 0, 69: 0}
    # the complete lists, written to slope_lists_le73.txt
    with open(os.path.join(HERE, "slope_lists_le73.txt"), "w", newline="\n") as fh:
        fh.write("# one-vertex sources: for every attained offset D <= 73, the odd parts d of the slopes with the least\n"
                 "# 2-exponent j (d:j means: d * 2^i is a slope of offset D, at some order of the base, iff i >= j).\n"
                 "# Complete at every order (THM-5129 (A), (B)); written by source_supply.py.\n")
        for D in sorted(prim):
            fh.write(f"{D}: " + " ".join(f"{dl}:{j}" for dl, j in sorted(prim[D].items())) + "\n")
    print(f"    complete lists for the {len(prim)} attained offsets <= {DMAX}: slope_lists_le73.txt "
          f"({sum(len(v) for v in prim.values())} odd parts in all)")
    # (4) two-step / near-join values below 3 * 75 = 225
    NMAX = 3 * (DMAX + 2)
    r = {}
    for N in range(3, NMAX, 2):
        c = 0
        for D in range(1, N // 3 + 1, 2):
            if is_slope(D, (N - D) // 2):
                c += 1
        r[N] = c
    zero = [N for N in r if r[N] == 0]
    print(f"(4) odd N in [3, {NMAX}) that are not D + 2U for any one-vertex source, at any order: {zero}")
    print("    number of representations r(N) = #{D : (D, (N - D)/2) a source}, N = 69, 71, ..., 131:")
    print("    " + " ".join(str(r[N]) for N in range(69, 133, 2)))
    lo = min(r[N] for N in range(129, NMAX, 2))
    print(f"    least r(N) on [129, {NMAX}): {lo}; mean {np.mean([r[N] for N in range(129, NMAX, 2)]):.2f}; "
          f"lattice points per line about N/6")
    ok &= zero == [7, 11, 13, 21, 47, 67]
    # (5) window sets of HYP-9453 (2)
    print("(5) h: number of slopes in [h, 3h/2] / number of integers there (attained h <= 73)")
    cells = []
    full_h = []
    for h in sorted(prim):
        hi = (3 * h) // 2
        tot = hi - h + 1
        c = sum(1 for u in range(h, hi + 1) if is_slope(h, u))
        cells.append(f"{h}:{c}/{tot}")
        if c == tot:
            full_h.append(h)
    line = "    "
    for c in cells:
        if len(line) + len(c) > 118:
            print(line)
            line = "    "
        line += c + " "
    print(line)
    print(f"    h with the whole window [h, 3h/2] in the slope set: {full_h}")
    ratio = [sum(1 for u in range(h, (3 * h) // 2 + 1) if is_slope(h, u)) / ((3 * h) // 2 - h + 1)
             for h in sorted(prim) if h >= 37]
    print(f"    mean covered fraction for 37 <= h <= 73: {np.mean(ratio):.3f} (min {min(ratio):.3f}, max {max(ratio):.3f})")
    # (6) strong classes per offset and their correlation with the covered share
    share = {}
    for h in sorted(prim):
        hi = (3 * h) // 2
        share[h] = sum(1 for u in range(h, hi + 1) if is_slope(h, u)) / (hi - h + 1)
    print("(6) h: number of strong isomorphism classes with H = h (orders 3..8; complete for h <= 73)")
    line = "    "
    for h in sorted(prim):
        c = f"{h}:{ncls.get(h, 0)} "
        if len(line) + len(c) > 118:
            print(line)
            line = "    "
        line += c
    print(line)
    for lo in (23, 15, 3):          # h = 1 (transitive, no strong class) is left out
        hs = [h for h in sorted(prim) if h >= lo]
        cc = np.corrcoef([ncls.get(h, 0) for h in hs], [share[h] for h in hs])[0, 1]
        print(f"    Pearson correlation of the covered share of [h, 3h/2] with the number of strong classes, "
              f"{lo} <= h <= 73 ({len(hs)} attained offsets): {cc:.2f}")
    poorest = sorted((share[h], h) for h in prim if h >= 23)[:6]
    print("    poorest offsets of [23, 73] (share, h, strong classes): "
          + ", ".join(f"({sh:.3f}, {h}, {ncls.get(h, 0)})" for sh, h in poorest))
    two = [h for h in sorted(prim) if ncls.get(h, 0) == 2]
    print(f"    offsets <= 73 with exactly two strong classes: {two}")
    ok &= two == [9, 17, 19, 35, 39, 47, 63]
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
