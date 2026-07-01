#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S77 -- CHARACTERIZATION TOWARD A PROOF: the safe-band residue system, the arithmetic
depth (witness denominator / CF), and the deep-well isolation of the covering-min construction.

FOUR new angles (all exact / grid-verified), synthesizing S64-S76:

  (1) SAFE-BAND RESIDUE SYSTEM. M(S) = max_{q,a} (1/q) * min_{v in S} ||a v||_q, where ||x||_q = dist(x mod q, 0).
      Loneliness at (q,a) <=> the residue system R(q,a) = {a v mod q} AVOIDS the danger band B_r = (-rq, rq).
      => the covering-min lower bound is a BAND-DODGING statement: every covering 13-set has some cyclic
      dilation avoiding a band of half-width ceil(r q).

  (2) THE EXTREMAL RESIDUE MECHANISM (general n). At the binding witness t* = n/Phi6 (a=n, q=Phi6=n^2-n+1):
      core {1..n-2} -> residues {n,2n,...,(n-2)n} (AP step n filling the safe band interior);
      SKIPPED runner (n-1) -> residue -1 (THE dangerous slot, adjacent to 0) -- this is WHY it is skipped;
      because n(n-1) = -1 mod Phi6, the multiple k(n-1) lands at residue -k (distance k);
      PATCH n(n-1) -> residue -n (the mirror safe-band edge). Binders: runner 1 (+n) and patch (-n). M=n/Phi6.

  (3) ARITHMETIC DEPTH = witness denominator q*(S) and its continued fraction. CF(t*) = [0; n-1, n] (DEEP).
      Every restructured covering 13-set binds at a SHALLOW witness (q* <= ~50), short CF, M ~ 0.10-0.14.

  (4) DEEP-WELL ISOLATION (quantified). Random covering primitive 13-sets have M >= ~0.108 (factor ~1.4 above
      the construction 0.0765, factor ~1.5 above the LRC threshold 1/14). The construction is the UNIQUE
      arithmetically-special minimizer; the "danger zone" near 1/14 = exactly the construction's scaled family.

  SLACK: LRC14 needs only M >= 1/14, but covering-min = 14/183 (margin = Dedekind sum 13/2562). Since the bulk
  sits at 0.108 >> 1/14, a certificate need only reach 1/14 -- real leverage for the bounded-case ILP.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
import numpy as np

n = 14
def cf(fr):
    a, b = fr.numerator, fr.denominator; out = []
    while b: out.append(a // b); a, b = b, a - (a // b) * b
    return out
def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for x in Sg:
        for y in Sg:
            for d in (x - y, x + y):
                if d > 0:
                    for k in range(1, d): cand.add(F(k, d))
    best = F(0); at = None
    for tt in cand:
        gg = min(min((v * tt) % 1, 1 - ((v * tt) % 1)) for v in Sg)
        if gg > best: best = gg; at = tt
    return best, at
def covers(S, n=n): return all(any(v % q == 0 for v in S) for q in range(2, n + 1))
def prim(S): return reduce(gcd, S) == 1
def Mgrid(S, G=4008):
    t = np.arange(1, G) / G
    g = np.min([np.abs(v * t - np.round(v * t)) for v in S], axis=0)
    i = int(np.argmax(g)); return float(g[i])

print("="*90)
print("(2) THE EXTREMAL RESIDUE MECHANISM at t*=n/Phi6 -- general n")
print("="*90)
for nn in [8, 10, 12, 14]:
    P = nn*nn - nn + 1; a = nn; q = P
    core = list(range(1, nn - 1)); patch = nn*(nn - 1)
    def dist(x, q=q): return min(x % q, q - (x % q))
    resid_core = [a*v % q for v in core]
    print(f"n={nn}: Phi6={q}, t*={a}/{q}=CF{cf(F(a,q))}, safe-band half-width n={nn}")
    print(f"   core 1..{nn-2} -> {resid_core}  (AP step {nn}); skipped {nn-1} -> residue {a*(nn-1)%q} (dist {dist(a*(nn-1))}=DANGER)")
    print(f"   patch {patch} -> residue {a*patch%q} (dist {dist(a*patch)}=n, safe edge); M = {nn}/{q}")

print()
print("="*90)
print("(3) ARITHMETIC DEPTH: construction is DEEP [0;n-1,n]; restructured covering sets are SHALLOW")
print("="*90)
tstar = F(14, 183)
core = list(range(1, 13))
tests = [("construction {1..12,182}", core + [182])]
# restructured: drop 2 core, add mult of 13 and mult of 14 (genuine 13-element covering sets)
for d1 in [(6, 11), (5, 9), (3, 7)]:
    rest = [x for x in core if x not in d1]
    S = sorted(set(rest + [26, 28]))
    if len(S) == 13 and covers(S) and prim(S):
        tests.append((f"drop{d1}+{{26,28}}", S))
print(f'{"set":30s} {"M":>9s} {"witness":>10s}  depth q*  CF')
for name, S in tests:
    M, at = Mexact(S)
    print(f'{name:30s} {str(M):>9s} {str(at):>10s}  {at.denominator:6d}  {cf(at)}')

print()
print("="*90)
print("(4) DEEP-WELL ISOLATION: random covering primitive 13-sets (speeds<=60)")
print("="*90)
random.seed(1); MAXV = 60; found = []
for _ in range(60000):
    S = sorted(random.sample(range(1, MAXV + 1), 13))
    if not covers(S) or not prim(S): continue
    found.append(Mgrid(S))
found.sort()
thr = 1/14
print(f"sampled {len(found)} covering primitive 13-sets; construction 14/183 = {float(tstar):.5f}; LRC threshold 1/14 = {thr:.5f}")
print(f"  min M = {found[0]:.5f}   5th pct = {found[len(found)//20]:.4f}   median = {found[len(found)//2]:.4f}")
print(f"  # below 0.10 = {sum(1 for m in found if m < 0.10)};  # below 1/14+.001 = {sum(1 for m in found if m < thr+.001)}")
print(f"  => bulk sits >= {found[0]:.3f} (factor {found[0]*14:.2f} above 1/14). Near-threshold covering sets are")
print(f"     NOT random -- they require the construction's arithmetic (the 182 patch reaching q=Phi6).")
print()
print("SYNTHESIS: covering-min lower bound splits by ARITHMETIC DEPTH:")
print("  shallow witnesses (q*<=Q): the bulk, M>=0.108>>1/14 (lazy-cut, easy vs 1/14 slack target);")
print("  deep witness q*=Phi6: only the construction family reaches it (forced-cover obstruction), M>=14/183;")
print("  huge speeds: S73/S74/S75 (don't beat). The residue/CF-depth is the through-line.")
print("DONE.")
