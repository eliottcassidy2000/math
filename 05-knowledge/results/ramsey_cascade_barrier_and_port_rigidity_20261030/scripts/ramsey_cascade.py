#!/usr/bin/env python3
"""ramsey_cascade.py -- the Greenwood-Gleason recursion P(k) = 2 - r + sum_i R(k with k_i - 1), its parity strictness,
Boza's mod-3 saturated-port criterion (arXiv:2603.10851, Thm 2.1) and the depth-t generalisation (K_t counts in the
nested saturated ports), evaluated mechanically on the DS1.18 (April 2026) table.  Also the exact requirements of the
owner's targets.  Conventions: R(k) for a sorted tuple k; entries with a 2 reduce to fewer colours; R(1, ...) = 1.
Inputs are CITED (Radziszowski DS1.18 sections 2.1 and 6.1; Boza 2026).  Standard library only."""
from math import factorial
from fractions import Fraction
import itertools, sys
# exact 2-colour values and best bounds (lo, hi) from DS1.18
two = {(3,3):(6,6),(3,4):(9,9),(3,5):(14,14),(3,6):(18,18),(3,7):(23,23),(3,8):(28,28),(3,9):(36,36),(3,10):(40,41),
       (4,4):(18,18),(4,5):(25,25),(4,6):(36,40),(4,7):(49,58),(4,8):(59,79),(5,5):(43,46),(5,6):(58,85),(5,7):(80,133),(6,6):(102,165)}
# multicolour known values / bounds (lo, hi)
multi = {(3,3,3):(17,17),(3,3,4):(30,30),(3,3,5):(45,57),(3,3,6):(61,91),(3,4,4):(55,77),(3,4,5):(89,157),(4,4,4):(128,229),
         (3,3,3,3):(51,62),(3,3,3,4):(97,149),(3,3,4,4):(174,450),(3,4,4,4):(381,1576),(4,4,4,4):(634,6301)}
def norm(k):
    k = tuple(sorted(x for x in k if x > 2))   # R(2, k, ...) = R(k, ...)
    return k
def bounds(k):
    """(lo, hi) for R(k) with k sorted, 2's removed; falls back to the recursion for the upper bound"""
    k = norm(k)
    if len(k) == 0: return (2, 2)   # only reached from tuples whose entries were all 2: R(2,...,2) = 2
    if len(k) == 1: return (k[0], k[0])
    if len(k) == 2:
        if k in two: return two[k]
        return (None, None)
    if k in multi: return multi[k]
    P = P_of(k)
    if P is None: return (None, None)
    U = [bounds(t)[1] for t in terms(k)]
    strict = (P % 2 == 0 and any(u % 2 == 0 for u in U))
    return (None, P - 1 if strict else P)
def terms(k):
    """the r terms R(k with k_i - 1), as sorted tuples (duplicates kept positionally)"""
    k = tuple(k); return [norm(k[:i] + (k[i]-1,) + k[i+1:]) for i in range(len(k))]
def P_of(k):
    k = norm(k); r = len(k); tot = 2 - r
    for t in terms(k):
        hi = bounds(t)[1]
        if hi is None: return None
        tot += hi
    return tot
def analyse(k):
    k = norm(k); r = len(k); T = terms(k); U = [bounds(t)[1] for t in T]
    if any(u is None for u in U): return None
    P = sum(U) + 2 - r
    out = {'k': k, 'terms': T, 'U': U, 'P': P, 'known': bounds(k)}
    out['parity_strict'] = (P % 2 == 0 and any(u % 2 == 0 for u in U))
    # Boza depth-1 and deeper cascades in each colour with k_i >= 4
    casc = []
    for i, ki in enumerate(k):
        if ki < 4: continue
        def level(j):
            kl = list(k); kl[i] = ki - j; return norm(kl)
        chain = []; note = None
        for j in range(1, ki - 1):
            tj = level(j); Uj = bounds(tj)[1]
            if Uj is None: note = 'missing input at %s' % (tj,); break
            if j >= 2:
                prev = level(j - 1); Up, Pp = bounds(prev)[1], P_of(prev)
                if Up is None or Pp is None or Up != Pp:
                    note = 'port %s unsaturated (U=%s, P=%s): deeper counts not forced' % (prev, Up, Pp); break
            chain.append(Uj - 1)
        n0 = P - 1; conds = []
        for t in range(3, len(chain) + 2):
            num = n0
            for j in range(t - 1): num *= chain[j]
            N_t = Fraction(num, factorial(t)); conds.append((t, chain[:t-1], N_t, N_t.denominator == 1))
        casc.append((i, chain, conds, note))
    out['cascades'] = casc
    # verdict
    strict_mod = any(any(not cc[3] for cc in c[2]) for c in casc)
    out['improves'] = out['parity_strict'] or strict_mod
    return out
def show(k):
    a = analyse(k)
    if a is None: print(k, 'inputs missing'); return
    print('R%s: known %s; P = %d from terms %s with U = %s%s' % (a['k'], a['known'], a['P'], a['terms'], a['U'], '; parity strict -> <= %d' % (a['P']-1) if a['parity_strict'] else ''))
    for (i, chain, conds, note) in a['cascades']:
        for (t, ports, N, isint) in conds:
            print('   colour %d (k=%d): nested port sizes %s -> forced number of colour-%d K_%d = %s: %s' % (i, a['k'][i], ports, i, t, N, 'integer' if isint else 'NOT an integer -> R <= P-1'))
        if note: print('   colour %d (k=%d): %s' % (i, a['k'][i], note))
    return a
print('=== Boza-style scan, DS1.18 inputs ===')
for k in sorted({tuple(sorted(t)) for t in itertools.product(range(3, 8), repeat=3)} | {tuple(sorted(t)) for t in itertools.product(range(3, 6), repeat=4)} | {(3,3,3,3,3),(3,3,3,3,4)}):
    show(k)
print()
print('=== What the owner targets require (recursion-level arithmetic) ===')
print("R(3,4,4) <= 73: the recursion gives 77 = 18 + 30 + 30 - 1 with all terms exact; no recursive route exists. Needs direct exclusion of (3,4,4;n) for n = 73..76.")
for x in (57, 56, 55, 54, 53, 52, 51):
    print('  if R(3,3,5) <= %d: P(3,3,6) = 2*18 + %d - 1 = %d%s' % (x, x, 35 + x, '' if (35+x) % 2 else ' (even; terms 18 even -> parity strict -> <= %d)' % (34 + x)))
print("R(3,3,6) <= 86 via recursion needs R(3,3,5) <= 52 with parity (P = 87 odd... ) -> check: R(3,3,5) <= 51 gives P = 86.")
for x, y in ((57,77),(56,77),(57,76),(51,73),(45,55)):
    P = 25 + x + y - 1; print('  R(3,4,5) <= %d from R(3,3,5) <= %d and R(3,4,4) <= %d%s' % (P, x, y, ' (parity strict -> <= %d)' % (P-1) if P % 2 == 0 and (x % 2 == 0 or y % 2 == 0) else ''))
print("R(3,4,5) <= 137 needs R(3,3,5) + R(3,4,4) <= 113 (current 134); e.g. 45 + 68, or 51 + 62.")
for y in (77, 76, 73, 70, 67, 66):
    P = 3 * y - 1; print('  R(4,4,4) <= %d from R(3,4,4) <= %d%s' % (P, y, ' (parity strict -> <= %d)' % (P-1) if P % 2 == 0 and y % 2 == 0 else ''))
print("R(4,4,4) <= 198 needs R(3,4,4) <= 66 (3*66 - 1 = 197).")
print()
print('=== (3,4,4;76) saturation ledger (exact consequences of saturation) ===')
n = 76; d1, d2, d3 = 17, 29, 29
print('degrees forced (17,29,29); colour-1 graph 17-regular triangle-free with', n*d1//2, 'edges; colour-2 and colour-3 graphs 29-regular K4-free with', n*d2//2, 'edges each')
print('N1(v) = unique (4,4;17) graph = Paley17 in colours 2/3: 68 colour-2 and 68 colour-3 edges, 8-regular each')
print('for r in N1(v): |N1(r) ∩ N2(v)| = |N1(r) ∩ N3(v)| = 8 (v has 8 colour-2 and 8 colour-3 neighbours in the Paley port N1(r)) -> e1(N1,N2) = e1(N1,N3) = 17*8 =', 17*8)
print('colour-1 degree sum over N2(v): 29*17 =', 29*17, '= 2 e1(N2) + 136 + e1(N2,N3); same for N3 -> e1(N2(v)) = e1(N3(v)) =: e*, 2e* + e1(N2,N3) =', 29*17 - 136)
print('   so e1(N2(v),N3(v)) is odd; with the in-port cap deg <= 8 (an (3,4)-set), e* <= 116 hence e1(N2,N3) >= %d' % (29*17 - 136 - 2*116))
print('mean colour-1 codegree |N1(w) ∩ N1(v)| over colour-2 edges vw: 136/29 =', Fraction(136, 29), '(non-integral: not constant)')
print('colour-2 triangles T2 = (1/3) sum_v e2(N2(v)); e2(N2(v)) <= 116; a mod-3 contradiction needs e2(N2(v)) forced to a constant c with 3 ∤ 76c -- NOT available: N2(v) is an unsaturated (3,3,4;29) port (30 < P(3,3,4) = 34)')
print('in-port degree caps in N2(v): (8, 8, 16) for colours (1, 2, 3), sum 28: slack 4 per vertex')
