"""Ternary coverage with the complete order-10 strong census (PREREGISTRATION.md, Q1-Q3).

Pieces (each an explicit strong tournament, so every product of piece values is realised by an order join, and H_3
is multiplicative over order joins, THM-4940 (1)):
  (i)   every strong class of order <= 9 (S86 census alpha_n.txt) and of order 10 (ternary_order10.c bitmaps);
  (ii)  every strong class with alpha_1 <= 80 (few_cycle_strong_80.txt, S86);
  (iii) strong one-apex tournaments (reduced words of length <= 19) and Session 16's two intercyclic rays.
The script forms the multiplicative closure of the piece values up to 1 + 3M, reports the threshold
m* = least alpha-sum m (m != 3, 10) with 1 + 3m outside the closure, and checks Q1 (m* >= 5000), Q2 (m* <= 20000) and
Q3 (every irreducible of M3 = {v = 1 mod 3} minus {10, 31} below 1 + 3m* is itself a piece value).
It also lists, lane by lane (ternary_lanes.py), the irreducibles of M3 below 1 + 3M that no piece attains.
This proves ATTAINMENT only; the only proved ternary gaps are 10 and 31 (THM-4940 (3), (4)).
usage: python ternary_coverage10.py CENSUS_DIR few_cycle_strong_80.txt T10_PREFIX M [pieces_out]
"""
import sys
import os
import glob
import importlib.util
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
s86 = os.path.join(here, "..", "packing_bases_repunits_20261108")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ta = load("ta", os.path.join(s86, "ternary_attainment.py"))
cp = ta.cp


def main():
    d, inv_path, t10, M = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    V = 1 + 3 * M
    vecs = ta.strong_census(d) | {tuple(v) for v in cp.load_inventory(inv_path)}
    apex = cp.one_apex_values(19, M)
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = {a for a in apex | rays if 1 <= a <= M}
    v_le9 = {1 + sum(c * 3 ** (k + 1) for k, c in enumerate(v)) for v in vecs}
    v_inter = {1 + 3 * a for a in inter}
    files = sorted(glob.glob(t10 + "_*.bin"))
    seen = None
    for f in files:
        b = np.fromfile(f, dtype=np.uint8)
        seen = b if seen is None else (seen | b)
    v_10 = set(int(x) for x in np.nonzero(seen)[0])
    print(f"order-10 bitmaps: {len(files)} slices, {len(v_10)} distinct strong H_3 values, max {max(v_10)}")
    pieces = v_le9 | v_inter | v_10
    bad = sorted(x for x in pieces if x % 3 != 1 or x in (10, 31))
    print(f"pieces: {len(v_le9)} values from orders <= 9 and alpha_1 <= 80, {len(v_inter)} intercyclic, "
          f"{len(v_10)} from order 10; {len(pieces)} distinct; values outside M3: {bad}")
    new10 = sorted(x for x in v_10 if x not in v_le9 and x not in v_inter)
    print(f"order-10 values not given by the S86 pieces: {len(new10)}")
    par, fac = ta.closure_with_witness([x for x in pieces if x <= V], V)
    unc = [m for m in range(0, M + 1) if par[1 + 3 * m] < 0]
    others = [m for m in unc if m not in (3, 10)]
    mstar = others[0] if others else None
    print(f"uncovered alpha-sums m <= {M}: {len(unc)}; m* = {mstar} (value {1 + 3 * mstar if mstar else None})")
    print(f"  first uncovered alpha-sums beyond 10: {others[:25]}")
    if mstar:
        print(f"Q1 (m* >= 5000): {'CONFIRMED' if mstar >= 5000 else 'REFUTED'}")
        print(f"Q2 (m* <= 20000): {'CONFIRMED' if mstar <= 20000 else 'REFUTED'}")
    for m in (11, 1800, 5000, (mstar or M) - 1):
        x = 1 + 3 * m
        f = []
        while x > 1:
            f.append(int(fac[x]))
            x = int(par[x])
        print(f"  witness m = {m}: {1 + 3 * m} = product of piece values {f}")
    # irreducibles of M3 below V by the reducibility sieve
    elems = np.array([v for v in range(4, V + 1) if v % 3 == 1 and v not in (10, 31)], dtype=np.int64)
    red = np.zeros(V + 1, dtype=bool)
    for x in elems:
        if x * x > V:
            break
        ys = elems[(elems >= x) & (elems <= V // x)]
        red[x * ys] = True
    irr = [int(v) for v in elems if not red[v]]
    if mstar:
        lim = 1 + 3 * mstar
        below = [v for v in irr if v < lim]
        miss_below = [v for v in below if v not in pieces]
        print(f"Q3: irreducibles of M3 below {lim}: {len(below)}; not a piece value: {len(miss_below)} "
              f"-> {'CONFIRMED' if not miss_below else 'REFUTED ' + str(miss_below[:10])}")
    miss = [v for v in irr if v not in pieces]
    by = {}
    for v in miss:
        by.setdefault(lane(v), []).append(v)
    print(f"irreducibles of M3 below {V}: {len(irr)}; not attained by any piece: {len(miss)}")
    for L in sorted(by):
        print(f"  lane {L}: {len(by[L])} missing; first {by[L][:12]}")
    if len(sys.argv) > 5:
        with open(sys.argv[5], "w") as fo:
            fo.write("\n".join(str(x) for x in sorted(pieces)) + "\n")


def factor(v):
    f = []
    d = 2
    while d * d <= v:
        while v % d == 0:
            f.append(d)
            v //= d
        d += 1
    if v > 1:
        f.append(v)
    return f


def lane(v):
    f = factor(v)
    r1 = [p for p in f if p % 3 == 1]
    r2 = [p for p in f if p % 3 == 2]
    if v in (40, 250, 310, 9610, 29791):
        return "sporadic"
    if len(f) == 1 and r1:
        return "P"
    if len(f) == 2 and len(r2) == 2:
        return "Q"
    if len(f) == 3 and sorted(r2) == [2, 5] and len(r1) == 1 and r1[0] != 31:
        return "10P"
    if len(f) == 2 and 31 in f and len(r1) == 2:
        return "31P"
    if len(f) == 3 and f.count(31) == 1 and len(r2) == 2:
        return "31Q"
    return "OTHER"


if __name__ == "__main__":
    main()
