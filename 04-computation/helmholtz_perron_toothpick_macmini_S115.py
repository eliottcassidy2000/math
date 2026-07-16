#!/usr/bin/env python3
"""helmholtz_perron_toothpick_macmini_S115.py -- mac-mini-2026-07-16-S115.
Owner brief: Helmholtz + arborescences beyond Hamiltonian paths; Kendall-Wei/Perron
self-composition invariant; toothpick/automata isomorphism hunting.

(A) THE HELMHOLTZ IDENTITY (exact): HodgeRank the unit arc flow of a tournament on K_n.
    Potential phi = d/(2n) (least squares); harmonic part = 0 (complete graph);
      GRADIENT ENERGY = sum_{u<v} (phi_u - phi_v)^2 = (n x)/(4 n^2) = x/(4n),
      CYCLE ENERGY   = C(n,2) - x/(4n).
    So THM-866's axis IS the Helmholtz gradient energy (x4n); the +8 walk is
    gradient-energy ascent; feedback (Kendall-Wei) lives on the curl side.
(B) KENDALL-WEI / PERRON census: lambda = spectral radius of the adjacency.
    (i) lambda = 0 <=> transitive (acyclic <=> nilpotent, one line) -- verify;
    (ii) corr(lambda, -x) over classes/labeled at n = 5, 6 (owner: ~0.93);
    (iii) lambda-spread within every x-level (strict refinement of the axis);
    (iv) corr(lambda, cycle-energy) -- the Helmholtz reading.
(C) TOOTHPICK probes: (i) A000568 = A002785 (mod 2) via the T <-> T^op involution
    (both even n >= 3); merged-metagraph counts V = (A000568+A002785)/2 =
    2, 3, 10, 34, 272, 3528, ... -- OEIS candidate string; (ii) v_2 ladder of A000568;
    (iii) generation-increment sequences for the growth-automaton framing.
"""
import sys, math
import numpy as np
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

print("(A) THE HELMHOLTZ IDENTITY -- exact check on random + exhaustive tournaments")
rng = np.random.default_rng(20260716)
ok = True
for trial in range(400):
    n = int(rng.integers(3, 10))
    A = np.zeros((n, n))
    for u, v in combinations(range(n), 2):
        if rng.random() < 0.5: A[u, v] = 1
        else: A[v, u] = 1
    s = A.sum(1)
    d = 2 * s - (n - 1)
    x = float((d ** 2).sum())
    phi = d / (2 * n)
    grad_e = sum((phi[u] - phi[v]) ** 2 for u, v in combinations(range(n), 2))
    ok &= abs(grad_e - x / (4 * n)) < 1e-9
print(f"   gradient energy == x/(4n) on 400 random tournaments n=3..9: {ok}")
print("   (proof: sum_(u<v) (d_u-d_v)^2 = n sum d^2 - (sum d)^2 = n x; /(2n)^2 = x/(4n))")
print("   cycle energy = C(n,2) - x/(4n); harmonic part = 0 (complete graph).")

print()
print("(B) KENDALL-WEI / PERRON")
for n in (5, 6):
    lams, xs, trans = [], [], 0
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    lam_by_x = {}
    for code in range(1 << m):
        A = np.zeros((n, n))
        for i, (u, v) in enumerate(pairs):
            if (code >> i) & 1: A[u, v] = 1
            else: A[v, u] = 1
        s = A.sum(1)
        d = 2 * s - (n - 1)
        x = int((d ** 2).sum())
        lam = max(abs(np.linalg.eigvals(A)))
        lams.append(lam); xs.append(x)
        lam_by_x.setdefault(x, []).append(lam)
        if lam < 1e-9:
            trans += 1
            ok_t = sorted(s) == list(range(n))
    lams = np.array(lams); xs = np.array(xs, dtype=float)
    corr = np.corrcoef(lams, -xs)[0, 1]
    spreads = {x: (min(v), max(v)) for x, v in sorted(lam_by_x.items())}
    nz = sum(1 for x, (lo, hi) in spreads.items() if hi - lo > 1e-9)
    print(f"   n={n}: labeled {1<<m}; lambda=0 count {trans} (= n! = {math.factorial(n)}: {trans == math.factorial(n)})"
          f"  corr(lambda, -x) = {corr:.4f}")
    print(f"        x-levels: {len(spreads)}; levels with lambda-spread > 0: {nz}"
          f"  (all nonterminal levels split: {nz >= len(spreads) - 2})")
    if n == 5:
        print("        per-level (x: lam-range):", {int(k): (round(v[0],3), round(v[1],3)) for k, v in list(spreads.items())})
    # Helmholtz reading: corr with cycle energy
    cyc = math.comb(n, 2) - xs / (4 * n)
    print(f"        corr(lambda, cycle-energy) = {np.corrcoef(lams, cyc)[0,1]:.4f}")

print()
print("(C) TOOTHPICK / AUTOMATA PROBES")
A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056, 11:903753248, 12:154108311168}
A002785 = {1:1, 2:1, 3:2, 4:2, 5:8, 6:12, 7:88, 8:176, 9:2752, 10:8784, 11:281504, 12:1414592}   # self-converse
par_ok = all((A000568[n] - A002785[n]) % 2 == 0 for n in A000568 if n in A002785)
print(f"   A000568 == A002785 (mod 2) for n <= 12: {par_ok}  [involution T <-> T^op: proved]")
merged = [ (A000568[n] + A002785[n]) // 2 for n in range(3, 13)]
print(f"   merged metagraph counts V(G_n/Z2) n=3..12: {merged}")
print(f"   OEIS candidate string: {','.join(map(str, merged))}")
v2 = [(n, (A000568[n] & -A000568[n]).bit_length() - 1) for n in range(3, 13)]
print(f"   v_2(A000568): {v2}")
inc = [A000568[n] - A000568[n-1] for n in range(4, 13)]
print(f"   generation increments (new mass per n): {inc[:6]}...")
print(f"   SC counts A002785 even for n >= 3 too: pairs within SC = quasi-fixed structure (THM-852 frame)")
print("\nDONE")
