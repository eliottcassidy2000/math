"""Packing counts in every base: H_b(T) = I(Omega(T), b) = sum_k alpha_k b^k, the number of b-coloured odd-cycle
packings (b = 2 is OCF: H_2 = H).  Exact gap tables for the values 1 + b*m with m <= MBUDGET, every base b.

Completeness argument (all steps elementary):
  * I(Omega(A => B), x) = I(Omega(A), x) I(Omega(B), x) (odd cycles never cross an order join), so the packing
    polynomial of any tournament is the product of those of its strong components.
  * If 1 + b*m = H_b(T) then alpha_1(T) <= m (all alpha_k >= 0 and alpha_1 is the coefficient of b^0 in m).
  * Moon's sharp cycle count (a strong k-tournament has at least k - l + 1 cycles of each length l, 3 <= l <= k;
    Moon, Canad. Math. Bull. 9 (1966) 297-301, Cor. 2.1, as cited in THM-343; MISTAKE-560) gives alpha_1 >= m(k) = sum over odd l in [3, k] of (k - l + 1): m(9) = 16, m(10) = 20.
    So with MBUDGET <= 19 every strong component of a realisation has order <= 9, and the census to order 9 is
    complete for this question.
Input (either):
  * a directory with alpha_<n>.txt for n = 3..9 written by alpha_census.c (complete for MBUDGET <= 19), or
  * a file written by few_cycle_strong.c M (every strong class with alpha_1 <= M, any order; complete for
    MBUDGET <= M by Moon's deletion lemma and sharp cycle count).
usage: python base_spectra.py CENSUS_DIR_or_FCS_FILE [MBUDGET=19] [BMAX=40]
"""
import sys
import os
from itertools import combinations_with_replacement


def moon_floor(k):
    return sum(k - l + 1 for l in range(3, k + 1, 2))


def load_fcs(path, mbudget):
    polys = {}
    for line in open(path):
        f = line.split()
        n = int(f[0]); a5 = [int(x) for x in f[1:6]]; H = int(f[6])
        # few_cycle_strong prints alpha_1..alpha_5 and H = sum_{k<=6} alpha_k 2^k; at order <= 18 at most six odd
        # cycles are pairwise disjoint, so alpha_6 is recovered exactly from H
        r = H - 1 - sum(c << (k + 1) for k, c in enumerate(a5))
        assert r >= 0 and r % 64 == 0 and n <= 18
        a = tuple(a5 + [r // 64])
        if a[0] <= mbudget:
            polys.setdefault(a, n)
    return polys


def load_strong_polys(d, mbudget):
    if os.path.isfile(d):
        return load_fcs(d, mbudget)
    polys = {}
    for n in range(3, 10):
        for line in open(os.path.join(d, f"alpha_{n}.txt")):
            f = line.split()
            strong = int(f[2])
            a = tuple(int(x) for x in f[4:8])
            if not strong or a[0] > mbudget:
                continue
            H = int(f[1])
            assert H == 1 + 2 * a[0] + 4 * a[1] + 8 * a[2] + 16 * a[3]
            polys.setdefault(a, n)
    return polys          # alpha-vector (a1..a4) -> least order


def mul(p, q):
    """multiply packing polynomials given as coefficient tuples (1, a1, a2, ...)"""
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i + j] += x * y
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    return tuple(r)


def closure_in_base(atoms, b, mbudget):
    """values 1 + b m (m <= mbudget) attained in base b: the multiplicative closure of the strong-atom values
    H_b(C) = I(Omega(C), b) (H_b is multiplicative over order joins).  Returns {m: witness atom list}."""
    top = 1 + b * mbudget
    vals = {}
    for a in atoms:
        v = 1 + sum(c * b ** (k + 1) for k, c in enumerate(a))
        if v <= top:
            vals.setdefault(v, a)
    reach = {1: ()}
    for v in sorted(vals):                       # unbounded knapsack over atoms, multiplicatively
        a = vals[v]
        for x in sorted(reach):
            y = x * v
            while y <= top:
                if y not in reach:
                    reach[y] = reach[x] + (a,)
                y *= v
                break
    # iterate to a fixed point (products of repeated atoms in any order)
    changed = True
    while changed:
        changed = False
        for v, a in vals.items():
            for x in list(reach):
                y = x * v
                if y <= top and y not in reach:
                    reach[y] = reach[x] + (a,)
                    changed = True
    return {(x - 1) // b: w for x, w in reach.items() if (x - 1) % b == 0}


def small_vectors(atoms, cap):
    """all packing polynomials (1, a1, a2, ...) of tournaments with alpha_1 <= cap (products of strong atoms)"""
    polys = [(1,) + a for a in atoms if a[0] <= cap]
    polys = [p[:max(i for i, c in enumerate(p) if c) + 1] for p in polys]
    reach = {(1,)}
    frontier = {(1,)}
    while frontier:
        new = set()
        for p in frontier:
            for q in polys:
                r = mul(p, q)
                if r[1] <= cap and r not in reach:
                    new.add(r)
        reach |= new
        frontier = new
    return reach


def decomposition_report(atoms, targets=(3, 10)):
    """for each target alpha-sum m and each base b <= m, every decomposition m = sum_k alpha_k b^(k-1) that respects
    alpha_k >= 1 => alpha_1 >= k, with its realisability"""
    cap = max(targets)
    V = small_vectors(atoms, cap)
    vecs = {tuple(list(p[1:]) + [0] * (6 - len(p[1:]))) for p in V}
    print(f"realisable packing vectors with alpha_1 <= {cap}: {len(V)}")
    for m in targets:
        print(f"  alpha-sum m = {m}:")
        for b in range(2, m + 2):
            decs = []
            def rec(k, rem, cur):
                if k == 6:
                    if rem == 0:
                        decs.append(tuple(cur))
                    return
                w = b ** k
                for c in range(rem // w + 1):
                    rec(k + 1, rem - c * w, cur + [c])
            rec(0, m, [])
            show = []
            for dvec in decs:
                if any(dvec[k] >= 1 and dvec[0] < k + 1 for k in range(1, 6)):
                    continue                       # alpha_{k+1} >= 1 forces alpha_1 >= k + 1
                tag = "REALISED" if dvec in vecs else "impossible"
                show.append(f"{tuple(x for x in dvec if True)[:3]} {tag}")
            print(f"    b = {b:2d}: " + "; ".join(show))


def main():
    d = sys.argv[1]
    mbudget = int(sys.argv[2]) if len(sys.argv) > 2 else 19
    bmax = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    if os.path.isdir(d):
        assert moon_floor(10) > mbudget, "census to order 9 is complete only for MBUDGET < m(10) = 20"
    strong = load_strong_polys(d, mbudget)
    atoms = sorted(strong)
    print(f"strong atoms (distinct packing vectors) with alpha_1 <= {mbudget}: {len(atoms)}")
    A = sorted({0} | {a[0] for a in atoms if all(c == 0 for c in a[1:])})
    print(f"A = alpha_1 of intercyclic tournaments (alpha_2 = 0), within budget: {A}")
    holes = [a for a in range(mbudget + 1) if a not in A]
    print(f"holes of A in [0, {mbudget}]: {holes}")
    print()
    print("decompositions of the two universal alpha-sums (vectors shown as (alpha_1, alpha_2, alpha_3)):")
    decomposition_report(atoms, (3, 10))
    print()
    print(f"gap table: for each base b, the m <= {mbudget} with 1 + b m NOT a b-coloured packing count")
    att = {}
    for b in range(2, bmax + 1):
        att[b] = closure_in_base(atoms, b, mbudget)
        g = [m for m in range(mbudget + 1) if m not in att[b]]
        print(f"  b = {b:2d}: missing m = {g}   (values 1 + b m = {[1 + b * m for m in g]})")
    print()
    print("rescue sets R(a) = {b <= BMAX : 1 + a b in Spec_b} for the holes a of A (for b > a only (a, 0) counts):")
    for a in holes:
        R = [b for b in range(2, bmax + 1) if a in att[b]]
        print(f"  a = {a:2d} (H = {1 + 2 * a:3d}): rescued in bases {R}")
        for b in R[:6]:
            print(f"      b = {b}: product of strong packing vectors {att[b][a]}")
    print()
    print("three-digit repunits 111_b = b^2 + b + 1 as b-coloured packing counts (b + 1 <= budget):")
    for b in range(2, min(mbudget, bmax + 1)):
        ok = (b + 1) in att[b]
        print(f"  b = {b:2d}: 111_b = {b * b + b + 1:5d}  {'attained' if ok else 'NOT attained'}  (b + 1 in A: {b + 1 in A})")


if __name__ == "__main__":
    main()
