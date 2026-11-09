"""Coverage of large alpha-sums by explicit products (PREREGISTRATION-2: P6-P8).

H_b is multiplicative over order joins, so every product of strong pieces is realised by a tournament (their order
join).  The pieces are:
  (i)   every strong class with alpha_1 <= 80, with its full packing vector (few_cycle_strong_80.txt);
  (ii)  strong one-apex tournaments: reduced words 1 u 0 of length <= LMAX, alpha = (a, 0), a = sum over i < j with
        w_i = 1, w_j = 0 of c(j - i), c(1) = 1, c(d) = 2^(d-2) (Session 16 / THM-4876 (6));
  (iii) Session 16's two intercyclic rays, alpha_1 = 3 * 2^s and 9 * 2^s - 2.
This proves ATTAINMENT only.  Gaps are proved only by the complete budget-80 tables (base_spectra_80.out), which this
script cross-checks on m <= 80.
usage: python coverage_products.py few_cycle_strong_80.txt
"""
import sys
import numpy as np


def one_apex_values(lmax, cap):
    vals = set()
    for m in range(2, lmax + 1):
        k = m - 2
        U = np.arange(1 << k, dtype=np.int64)
        if k:
            bits = ((U[:, None] >> np.arange(k - 1, -1, -1)) & 1).astype(np.int64)
        else:
            bits = np.zeros((1, 0), np.int64)
        W = np.concatenate([np.ones((U.size, 1), np.int64), bits, np.zeros((U.size, 1), np.int64)], axis=1)
        val = np.zeros(U.size, dtype=np.int64)
        for i in range(m):
            for j in range(i + 1, m):
                d = j - i
                c = 1 if d == 1 else 1 << (d - 2)
                val += W[:, i] * (1 - W[:, j]) * c
        vals.update(int(v) for v in np.unique(val[val <= cap]))
    return vals


def load_inventory(path):
    vecs = set()
    for line in open(path):
        f = line.split()
        n = int(f[0])
        a5 = [int(x) for x in f[1:6]]
        H = int(f[6])
        r = H - 1 - sum(c << (k + 1) for k, c in enumerate(a5))
        assert r >= 0 and r % 64 == 0 and n <= 18
        vecs.add(tuple(a5 + [r // 64]))
    return vecs


def closure(atom_values, V):
    reach = np.zeros(V + 1, dtype=bool)
    reach[1] = True
    for t in sorted(set(atom_values)):
        if t > V or t < 2:
            continue
        xs = np.arange(1, V // t + 1)
        while True:
            src = reach[xs]
            tgt = xs * t
            new = src & ~reach[tgt]
            if not new.any():
                break
            reach[tgt[new]] = True
    return reach


def main():
    inv = load_inventory(sys.argv[1])
    plan = {2: 100000, 3: 100000, 4: 20000, 5: 20000, 6: 20000, 7: 20000, 8: 20000}
    amax = max(plan.values())
    apex = one_apex_values(21, amax)               # reduced length 21 gives every one-apex value below 2^18
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = {a for a in apex | rays if 1 <= a <= amax}
    print(f"pieces: {len(inv)} strong packing vectors with alpha_1 <= 80; {len(inter)} intercyclic alpha_1 values <= {amax}")
    exact_gaps = {2: [3, 10], 3: [3, 10], 4: [3, 14, 27], 5: [3, 10, 27], 6: [3, 10, 14, 21, 46], 7: [3, 10, 14, 29],
                  8: [3, 35, 52]}
    for b, M in plan.items():
        V = 1 + b * M
        vals = [1 + sum(c * b ** (k + 1) for k, c in enumerate(v)) for v in inv]
        vals += [1 + b * a for a in inter]
        reach = closure([v for v in vals if v <= V], V)
        low = [m for m in range(0, 81) if not reach[1 + b * m]]
        unc = [m for m in range(81, M + 1) if not reach[1 + b * m]]
        print(f"b = {b}: uncovered alpha-sums m <= 80: {low} (exact gap table: {exact_gaps[b]}, "
              f"{'agrees' if low == exact_gaps[b] else 'DISAGREES'}); "
              f"uncovered in (80, {M}]: {len(unc)} {unc[:20]}")


if __name__ == "__main__":
    main()
