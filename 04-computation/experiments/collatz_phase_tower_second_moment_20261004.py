#!/usr/bin/env python3
"""
E5o (opus, 2026-10-04): the exact second moment of the random-start Pascal tower (Fourier note 4k).

Model: the deepest-level digit string is a uniform residue R mod 2^M (the start x of the tower x 3^-n); all shallower
levels follow by the Pascal identity theta_{n-m,d} = sum_i C(m,i) theta_{n,d-i}, i.e. the exact phases of the
Collatz tower with a random 2-adic start.  Then
    f_n(0) = sum_paths 2^-|a| e( sum_i theta_{i, D_i} ),   D_i = a_i + ... + a_n,
and, because theta_{i,D} = 3^(n-i) theta_{n,D} = 3^(n-i) R / 2^D (mod 1), the path phase is e(R Phi(a)) with the REAL
dyadic value Phi(a) = sum_i 3^(n-i) 2^-D_i (mod 1).  Averaging over R uniform mod 2^M (M >= max depth):
    E_R |f_n(0)|^2 = sum_{a,b} 2^-|a|-|b| 1[Phi(a) = Phi(b) mod 1]  =  sum_{groups} ( sum_{a in group} 2^-|a| )^2  >=  3^-n,
with equality iff no two distinct paths collide.  (For R restricted to odd residues the pairs with Phi(a) - Phi(b) = 1/2
contribute -1 instead of 0.)  We compute the group sums exactly (valuations a_i <= A) and compare with (a) the i.i.d.
model's exact 3^-n, (b) a Monte-Carlo over random starts R (all residues) and over random odd starts, and (c) the
collision mass by n.
Usage: python collatz_phase_tower_second_moment_20261004.py [NMAX=6] [A=12]
"""
import sys, math, time
from itertools import product
from collections import defaultdict
import numpy as np

def exact_second_moment(n, A):
    """returns (sum of squared group weights, diagonal 3^-n-type sum, number of groups, collision mass) at level n."""
    Dmax = n * A
    groups = defaultdict(float)
    diag = 0.0
    scale = 1 << Dmax
    for a in product(range(1, A + 1), repeat=n):          # a = (a_1, ..., a_n)
        w = 2.0 ** (-sum(a))
        # Phi(a) * 2^Dmax as an integer mod 2^Dmax: sum_i 3^(n-i) 2^(Dmax - D_i), D_i = a_i + ... + a_n
        key = 0; D = 0
        for i in range(n, 0, -1):                          # i = n down to 1
            D += a[i - 1]
            key += 3 ** (n - i) * (1 << (Dmax - D))
        key %= scale
        groups[key] += w
        diag += w * w
    tot = sum(v * v for v in groups.values())
    return tot, diag, len(groups), tot - diag

def mc_second_moment(n, A, M, S, odd_only, seed=0):
    """Monte Carlo over random starts: f_n(0) computed directly from the exact phases of the random start."""
    rng = np.random.default_rng(seed)
    vals = []
    paths = list(product(range(1, A + 1), repeat=n))
    ws = np.array([2.0 ** (-sum(a)) for a in paths])
    # D_i for each path and the real dyadic Phi(a) as a float (M <= 50 bits keeps float exact enough)
    Phi = []
    for a in paths:
        D = 0; ph = 0.0
        for i in range(n, 0, -1):
            D += a[i - 1]
            ph += 3 ** (n - i) * 2.0 ** (-D)
        Phi.append(ph % 1.0)
    Phi = np.array(Phi)
    for s in range(S):
        R = int(rng.integers(0, 1 << M))
        if odd_only: R |= 1
        f = np.sum(ws * np.exp(2j * np.pi * ((R * Phi) % 1.0)))
        vals.append(abs(f) ** 2)
    return float(np.mean(vals)), float(np.std(vals) / math.sqrt(S))

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"exact second moment of the random-start Pascal tower, valuations <= {A}")
    P("n | 3^n E_R|f_n(0)|^2 (exact group sum) | 3^n x diagonal (iid value, -> 1 minus truncation) | #groups / #paths | 3^n x collision mass | MC all residues | MC odd residues")
    t0 = time.time()
    for n in range(1, NMAX + 1):
        tot, diag, ng, coll = exact_second_moment(n, A)
        npaths = A ** n
        mc_all = mc_second_moment(n, min(A, 8), 40, 400, False) if n <= 5 else (float("nan"), 0)
        mc_odd = mc_second_moment(n, min(A, 8), 40, 400, True) if n <= 5 else (float("nan"), 0)
        P(f"{n} | {3**n * tot:.6f} | {3**n * diag:.6f} | {ng} / {npaths} | {3**n * coll:.6f} | {3**n * mc_all[0]:.3f} +- {3**n * mc_all[1]:.3f} (A<=8) | {3**n * mc_odd[0]:.3f} +- {3**n * mc_odd[1]:.3f}   [{time.time()-t0:.0f}s]")
    P("reading: the exact group sum exceeds the diagonal by the collision mass (>= 0): the random-start Collatz tower has a LARGER mean square than the i.i.d. model, never smaller;")
    P("the Monte-Carlo over all residues should match the exact group sum (truncation A <= 8 there), the odd-only one may sit lower by the half-collision pairs.")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
