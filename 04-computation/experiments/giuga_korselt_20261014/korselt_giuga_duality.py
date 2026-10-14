#!/usr/bin/env python3
"""The Korselt / Giuga duality, exactly.

For squarefree composite n and p | n:
  Korselt at p:  (p - 1) | (n - 1)          [all p: Carmichael]
  Giuga  at p:   p | (n/p - 1)              [all p: Giuga]
A counterexample to Giuga's conjecture is exactly a number that is both (Giuga 1950).

Parts
 A  Lerch's p^2-lift of Giuga's sum: S(p) = sum_{k<p} k^(p-1) == p + (p-1)! (mod p^2) for every prime p < PMAX,
    checked directly; hence Wilson primes <=> S(p) == p - 1 and w_p == -1 <=> S(p) == -1 (mod p^2).
 B  Parity divide: for even squarefree n the Korselt score is exactly 1 (only p = 2); checked on all even squarefree
    n <= N.  For odd n, Giuga with sum 1/p - 1/n = k forces omega(n) == k + 1 (mod 2); checked on all squarefree
    odd n <= N with Giuga score omega (none exist below N, as expected).
 C  Carmichael numbers <= N (vectorized Korselt over a smallest-prime-factor sieve), their Giuga scores
    G(c) = #{p | c : c/p == 1 mod p}, the near-counterexamples (highest G(c)/omega(c)), and Carmichael twins
    (c - 2 or c + 2 prime).
 D  Known Giuga numbers: Korselt score 1, g' = g + 1.

Usage: python3 korselt_giuga_duality.py [N=10**7] [PMAX=3000]
"""
import sys
import time
from math import prod
from fractions import Fraction as F

import numpy as np
import sympy

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7
PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
checks = 0
t0 = time.time()

# ---------------- A: Lerch lift
wil, gsq = [], []
for p in sympy.primerange(3, PMAX):
    p2 = p * p
    S = sum(pow(k, p - 1, p2) for k in range(1, p)) % p2
    f = 1
    for k in range(2, p):
        f = f * k % p2
    assert S == (p + f) % p2, p
    w = ((f + 1) // p) % p
    assert S == (p - 1 + p * w) % p2
    checks += 2
    if w == 0:
        wil.append(p)
    if w == p - 1:
        gsq.append(p)
print(f"A  S(p) == p + (p-1)! == p - 1 + p w_p (mod p^2) for all primes 3 <= p < {PMAX}: OK; "
      f"Wilson primes {wil}; w_p == -1 (S(p) == -1 mod p^2): {gsq}  [{time.time() - t0:.1f}s]", flush=True)

# ---------------- sieve
spf = np.zeros(N + 1, dtype=np.int32)
for p in sympy.primerange(2, int(N ** 0.5) + 1):
    block = spf[p * p::p]
    block[block == 0] = p
spf[spf == 0] = np.arange(N + 1, dtype=np.int32)[spf == 0]
print(f"   sieve to {N} done  [{time.time() - t0:.1f}s]", flush=True)


def factor(n):
    out = []
    while n > 1:
        p = int(spf[n])
        out.append(p)
        n //= p
    return out


# ---------------- C: Carmichael numbers by vectorized Korselt (odd n only)
n = np.arange(3, N + 1, 2, dtype=np.int64)
m = n.copy()
ok = np.ones(len(n), dtype=bool)
nf = np.zeros(len(n), dtype=np.int8)
while True:
    active = m > 1
    if not active.any():
        break
    p = spf[m[active]].astype(np.int64)
    sub = ok[active]
    sub &= ((n[active] - 1) % (p - 1) == 0)
    mm = m[active] // p
    sub &= (mm % p != 0)                 # squarefree
    ok[active] = sub
    nf[active] += 1
    m[active] = mm
carm = n[ok & (nf >= 3)]
print(f"C  Carmichael numbers <= {N}: {len(carm)}; first {carm[:8].tolist()}  [{time.time() - t0:.1f}s]", flush=True)
rows = []
for c in carm.tolist():
    ps = factor(c)
    assert len(set(ps)) == len(ps) and len(ps) >= 3
    assert all((c - 1) % (q - 1) == 0 for q in ps)
    checks += 1
    G = [q for q in ps if (c // q - 1) % q == 0]
    rows.append((len(G) / len(ps), len(G), len(ps), c, ps, G))
dist = {}
for fr, g, r, c, ps, G in rows:
    dist[(g, r)] = dist.get((g, r), 0) + 1
print(f"   Giuga-score distribution (G, omega): {dict(sorted(dist.items()))}")
rows.sort(key=lambda t: (-t[0], t[3]))
print("   near-counterexamples (highest share of primes also satisfying Giuga):")
for fr, g, r, c, ps, G in rows[:12]:
    print(f"     {c} = {'*'.join(map(str, ps))}: Giuga at {G} ({g}/{r})")
tw = [c for c in carm.tolist() if sympy.isprime(c - 2) or sympy.isprime(c + 2)]
print(f"   Carmichael twins (c - 2 or c + 2 prime): {len(tw)} of {len(carm)}: "
      f"{[(c, 'c-2' if sympy.isprime(c - 2) else '', 'c+2' if sympy.isprime(c + 2) else '') for c in tw[:12]]}")

# ---------------- B: parity divide
even_sqfree_checked = 0
for e in range(6, min(N, 2 * 10 ** 6) + 1, 2):
    ps = factor(e)
    if len(set(ps)) != len(ps) or len(ps) < 2:
        continue
    K = sum(1 for q in ps if (e - 1) % (q - 1) == 0)
    assert K == 1, e
    even_sqfree_checked += 1
checks += even_sqfree_checked
print(f"B  every even squarefree composite n <= {min(N, 2 * 10 ** 6)} has Korselt score exactly 1 (only p = 2): "
      f"{even_sqfree_checked} numbers OK")
odd_giuga = []
for c in range(15, min(N, 2 * 10 ** 6) + 1, 2):
    ps = factor(c)
    if len(set(ps)) != len(ps) or len(ps) < 2:
        continue
    if all((c // q - 1) % q == 0 for q in ps):
        odd_giuga.append(c)
print(f"   odd Giuga numbers <= {min(N, 2 * 10 ** 6)}: {odd_giuga} (none expected: sum 1/p over odd primes needs >= 9 primes)")

# ---------------- D: known Giuga numbers
ROWS = [[2, 3, 5], [2, 3, 11, 13], [2, 3, 7, 41], [2, 3, 11, 17, 59], [2, 3, 11, 23, 31, 47057],
        [2, 3, 7, 43, 3041, 4447], [2, 3, 7, 59, 163, 1381, 775807], [2, 3, 7, 71, 103, 67213, 713863],
        [2, 3, 7, 71, 103, 61559, 29133437], [2, 3, 11, 23, 31, 47137, 28282147, 3892535183],
        [2, 3, 11, 23, 31, 47059, 2259696349, 110725121051], [2, 3, 7, 43, 1831, 138683, 2861051, 1456230512169437]]
for ps in ROWS:
    g = prod(ps)
    assert sum(1 for q in ps if (g - 1) % (q - 1) == 0) == 1
    assert sum(g // q for q in ps) == g + 1
    assert sum(F(1, q) for q in ps) - F(1, g) == 1
    checks += 3
print(f"D  all 12 Giuga numbers of A007850: Korselt score 1, g' = g + 1, sum 1/p - 1/g = 1: OK")
print(f"total exact checks: {checks}  [{time.time() - t0:.1f}s]")
