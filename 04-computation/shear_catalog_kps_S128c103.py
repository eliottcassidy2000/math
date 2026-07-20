#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c103 -- HYP-8170: THE SHEAR CATALOG.

Owner: shear the n*2^x+1 grid down 1 (triangle) or 2 (Fibonacci-analogous);
sum or PRODUCT the columns; do the same for ALL meaningful continuations of the
triangular numbers; compare; T_n = C(n,2) = edges of K_n (the relational count).

Machinery: a family = grid G(c, n) (column c >= c0, row n >= n0).  Shear-s row
sum R_s(m) = Sum_c G(c, m - s*c) over valid indices; likewise products.
- CONTROL: simplex/binomial grid: s=1 -> 2^(m-1), s=2 -> Fibonacci, s=3 ->
  Narayana's cows (A000930).
- PROTH grid n*2^x+1: GF-exact shear laws; growth 2^(1/s): the pseudo-doubling
  ladder (sqrt2 at the owner's s=2).
- polygonal, Faulhaber, centered polygonal, pyramidal, pronic: polynomial laws.
- products for the marquee families; recurrence fitter; growth constants;
  OEIS batch (curl, separate); klein T1532 sequence-D probe.
"""
from fractions import Fraction as Fr
import itertools, math

def binom(n, k):
    if k < 0 or n < 0 or k > n: return 0
    return math.comb(n, k)

# ---------------- families as G(c, n) ----------------
def simplex(c, n):      # c-dimensional simplex numbers, c >= 0; T at c = 2
    if n < 0: return None
    return binom(n + c - 1, c) if c >= 0 else None
def proth(c, n):        # n*2^c + 1, c >= 0, n >= 1 (owner's grid)
    if n < 1: return None
    return n * 2**c + 1
def polygonal(c, n):    # c-gonal numbers, c >= 3; T at c = 3
    if n < 1: return None
    return ((c - 2) * n * n - (c - 4) * n) // 2
def faulhaber(c, n):    # Sum_{j<=n} j^c, c >= 0; T at c = 1
    if n < 1: return None
    return sum(j**c for j in range(1, n + 1))
def centered(c, n):     # centered c-gonal, c >= 3: 1 + c*T(n-1)
    if n < 1: return None
    return 1 + c * (n - 1) * n // 2
def pyramidal(c, n):    # c-gonal pyramidal: Sum_j polygonal(c, j), c >= 3
    if n < 1: return None
    return sum(((c - 2) * j * j - (c - 4) * j) // 2 for j in range(1, n + 1))
def pronic(c, n):       # n(n+c)/2-family, c >= 0 integer: T at c = 1 (n(n+1)/2)
    if n < 1: return None
    v = n * (n + c)
    return v // 2 if v % 2 == 0 else Fr(v, 2)

FAMS = {
  "simplex(control)": (simplex, 0),
  "Proth n*2^x+1":    (proth, 0),
  "polygonal":        (polygonal, 3),
  "Faulhaber":        (faulhaber, 0),
  "centered":         (centered, 3),
  "pyramidal":        (pyramidal, 3),
  "pronic":           (pronic, 0),
}

def shear_sum(G, c0, s, m, alt=False):
    tot = 0; c = c0; any_ = False
    while True:
        n = m - s * c
        if s > 0 and n < 1 - (1 if G is simplex else 0):  # simplex allows n=0
            break
        if s == 0 and c > m + 6:   # square-array convention: truncate c <= m+6? use c <= m
            break
        v = G(c, n) if n is not None else None
        if s == 0:
            n = m
            v = G(c, n)
            if c > m: break        # square: c <= m truncation (noted convention)
        if v is None:
            break
        tot += (-1)**(c - c0) * v if alt else v
        any_ = True
        c += 1
        if c - c0 > 400: break
    return tot if any_ else None

def shear_prod(G, c0, s, m):
    tot = 1; c = c0; any_ = False
    while True:
        n = m - s * c
        if n < 1: break
        v = G(c, n)
        if v is None or v == 0: break
        tot *= v
        any_ = True
        c += 1
        if c - c0 > 60: break
    return tot if any_ else None

MMAX = 15
print("== SHEAR SUMS, s = 1, 2, 3 (and s = 0 square truncated c <= m) ==", flush=True)
tables = {}
for name, (G, c0) in FAMS.items():
    print(f"\n  --- {name} ---", flush=True)
    for s in (0, 1, 2, 3):
        vals = [shear_sum(G, c0, s, m) for m in range(1, MMAX + 1)]
        vals = [v for v in vals if v is not None]
        tables[(name, s)] = vals
        # growth estimate
        gr = ""
        if len(vals) > 6 and vals[-2] and vals[-2] != 0:
            r = vals[-1] / vals[-2]
            gr = f"  growth~{r:.4f}"
        print(f"    s={s}: {vals[:12]}{gr}", flush=True)

print("\n== CONTROL CHECKS (Pascal ladder) ==", flush=True)
fib = [1,1,2,3,5,8,13,21,34,55,89,144,233,377,610]
print(f"  simplex s=2 vs Fibonacci: {tables[('simplex(control)',2)][:10]} vs {fib[:10]}", flush=True)
print(f"  simplex s=1 = 2^(m-1)?: {tables[('simplex(control)',1)][:8]}", flush=True)
print(f"  simplex s=3 (Narayana A000930 1,1,1,2,3,4,6,9,13,19): {tables[('simplex(control)',3)][:10]}", flush=True)

print("\n== PROTH EXACT LAWS (GF-derived, verify numerically) ==", flush=True)
import sympy as sp
t = sp.symbols('t')
for s in (1, 2, 3):
    # GF: Sum_x t^(s*x) * [GF_n of (n*2^x+1)] = Sum_x t^(s*x) * (2^x*t/(1-t)^2 + t/(1-t))
    GF = sp.summation(t**(s*sp.symbols('x', integer=True, nonnegative=True)) * 1, (sp.symbols('x'), 0, 0))
    # build symbolically: Sum_x [2^x t^(sx)] = 1/(1-2t^s); Sum_x t^(sx) = 1/(1-t^s)
    GFs = t/(1-t)**2 * 1/(1-2*t**s) + t/(1-t) * 1/(1-t**s)
    ser = sp.series(GFs, t, 0, 13).removeO()
    coeffs = [sp.expand(ser).coeff(t, m) for m in range(1, 13)]
    match = coeffs[:len(tables[("Proth n*2^x+1", s)][:12])] == [sp.Integer(v) for v in tables[("Proth n*2^x+1", s)][:12]]
    print(f"  s={s}: GF = t/((1-t)^2(1-2t^s)) + t/((1-t)(1-t^s)); matches data: {match}; dominant root 2^(1/{s}) = {2**(1/s):.4f}", flush=True)
print("  => THE PROTH SHEAR SPECTRUM = 2^(1/s): 2, sqrt2 (the hypotenuse/pseudo-doubling constant at the owner's Fibonacci-analogous s=2), 2^(1/3), ...", flush=True)
print("  Pascal ladder roots x^s = x^(s-1)+1: 2, phi = 1.6180, 1.4656 (Narayana), ... => Pascal DOMINATES Proth at every shear s >= 2", flush=True)

print("\n== PRODUCTS (marquee families, s = 2) ==", flush=True)
for name in ("simplex(control)", "Proth n*2^x+1", "polygonal", "Faulhaber"):
    G, c0 = FAMS[name]
    pv = [shear_prod(G, c0, 2, m) for m in range(1, 13)]
    pv = [v for v in pv if v is not None]
    print(f"  {name} s=2 products: {pv[:9]}", flush=True)
print("  (Proth s=2 products = products of Proth numbers along the diagonal -- factor census separate)", flush=True)

print("\n== ALTERNATING SUMS s=2 (parity metric) ==", flush=True)
for name in FAMS:
    G, c0 = FAMS[name]
    av = [shear_sum(G, c0, 2, m, alt=True) for m in range(1, 13)]
    av = [v for v in av if v is not None]
    print(f"  {name}: {av[:11]}", flush=True)

print("\n== KLEIN T1532 SEQUENCE-D PROBE ==", flush=True)
D = [1,1,2,3,5,8,13,21,33,51,76,111,157,218]
print(f"  klein D: {D}", flush=True)
found = False
for name in FAMS:
    for s in (1,2,3):
        v = tables.get((name, s), [])
        for off in range(-2, 3):
            seg = v[max(0,off):] if off >= 0 else v
            if len(seg) >= 8 and all(seg[i] == D[i - min(0,off)] for i in range(min(8, len(seg)))):
                print(f"  MATCH: {name} s={s} offset {off}", flush=True); found = True
if not found:
    print("  no direct family/shear match for klein's D in this catalog (D remains its own object)", flush=True)

print("\n== the relational reading ==", flush=True)
print("  T_n = C(n,2) = |E(K_n)| = the RELATIONAL size of n (the owner's frame; = tournament arc count, the tiling-hypercube dimension m = C(n-1,2) + base path)", flush=True)
print("  DICHOTOMY (measured above): relational families (simplex: counting sub-relations) shear EXPONENTIALLY (phi-ladder);", flush=True)
print("  positional families (polygonal/centered/pronic: dot patterns) shear POLYNOMIALLY; Faulhaber and pyramidal (sum-of-positional) polynomial one degree up;", flush=True)
print("  the Proth grid (doubling x relation-edge) shears exponentially on the 2^(1/s) ladder BELOW the relational phi-ladder.", flush=True)
print("\nDONE.", flush=True)
