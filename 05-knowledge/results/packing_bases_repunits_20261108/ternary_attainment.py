"""Constructive attainment in bases 2 and 3 from explicit products (post-audit addendum to THM-4940 / HYP-9530 (A)).

Pieces (each an explicit strong tournament, so every product is realised by an order join):
  (a) every strong class of order <= 9 (complete census from kind-pasteur's classgen.c + alpha_census.c);
  (b) every strong class with alpha_1 <= 80 (few_cycle_strong_80.txt);
  (c) strong one-apex tournaments (reduced words of length <= 19) and Session 16's two intercyclic rays.
For base b the script forms the multiplicative closure of the piece values H_b up to 1 + b*M and reports the uncovered
alpha-sums.  Attainment is witnessed by the closure (a factorisation is printed for spot checks).  Gaps are NOT
proved here: the only proved gaps are THM-4940 (3) and (4) (alpha-sums 3 and 10 in bases 2 and 3) and the complete
tables of THM-4940 (6).
This range was found by an exploratory calibration after PREREGISTRATION-2's P6-P8 failed; it is a constructive
verification, not a pre-registered test.
usage: python ternary_attainment.py CENSUS_DIR few_cycle_strong_80.txt
"""
import sys
import importlib.util
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("cp", os.path.join(here, "coverage_products.py"))
cp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cp)


def strong_census(d, maxn=9):
    vecs = set()
    for n in range(3, maxn + 1):
        for line in open(os.path.join(d, f"alpha_{n}.txt")):
            f = line.split()
            if f[2] == '1':
                vecs.add(tuple(int(x) for x in f[4:8]))
    return vecs


def closure_with_witness(atom_values, V):
    """multiplicative closure; par[x] = (y, t) with x = y * t"""
    par = np.full(V + 1, -1, dtype=np.int64)
    fac = np.zeros(V + 1, dtype=np.int64)
    par[1] = 0
    for t in sorted(set(atom_values)):
        if t > V or t < 2:
            continue
        xs = np.arange(1, V // t + 1)
        while True:
            src = par[xs] >= 0
            tgt = xs * t
            new = src & (par[tgt] < 0)
            if not new.any():
                break
            par[tgt[new]] = xs[new]
            fac[tgt[new]] = t
    return par, fac


def main():
    d, inv_path = sys.argv[1], sys.argv[2]
    vecs = strong_census(d) | {v[:4] for v in cp.load_inventory(inv_path)}
    apex = cp.one_apex_values(19, 70000)
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = {a for a in apex | rays if 1 <= a <= 70000}
    print(f"pieces: {len(vecs)} strong packing vectors (orders <= 9, plus alpha_1 <= 80), {len(inter)} intercyclic values")
    for b, M in ((2, 10000), (3, 10000)):
        V = 1 + b * M
        vals = [1 + sum(c * b ** (k + 1) for k, c in enumerate(v)) for v in vecs] + [1 + b * a for a in inter]
        par, fac = closure_with_witness([x for x in vals if x <= V], V)
        unc = [m for m in range(0, M + 1) if par[1 + b * m] < 0]
        first_other = next((m for m in unc if m not in (3, 10)), None)
        print(f"base {b}: uncovered alpha-sums m <= {M}: {len(unc)}; the only ones below {first_other} are "
              f"{[m for m in unc if first_other is None or m < first_other]}")
        print(f"  so every value 1 + {b}m with m < {first_other}, m not 3 or 10, is a base-{b} packing count "
              f"(values up to {1 + b * (first_other - 1)})")
        for m in (11, 100, 1000, first_other - 1):
            x = 1 + b * m
            f = []
            while x > 1:
                f.append(int(fac[x]))
                x = int(par[x])
            print(f"  witness m = {m}: 1 + {b}*{m} = {1 + b * m} = product of piece values {f}")


if __name__ == "__main__":
    main()
