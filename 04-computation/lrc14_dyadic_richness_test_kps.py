import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

# ============================================================================
# DYADIC RICHNESS HYPOTHESIS for LRC(14) sector route.
#
# p_0(E) = meas(S7(E)) = meas{ x in [0,1): orbit {frac(e x): e in E u {0}}
#                              hits all 7 sectors [j/7,(j+1)/7) }.
# E = speed set (positive ints). We test whether p_0 is governed by the
# "dyadic richness" of E: how much 2-adic doubling-chain structure E carries.
#
# Glaisher: each e = 2^a * b (b odd). The orbit decomposes as a union over
# odd b of dyadic doubling-orbits {frac(2^a b x): a}. The doubling map z->2z
# acts with order 3 on sectors mod 7, splitting inner sectors into QR={1,2,4}
# and NQR={3,5,6}. Richer dyadic chains => more sectors hit per "free" speed.
# ============================================================================

def meas_S7(E):
    """EXACT meas{x in [0,1): {frac(e x): e in E u {0}} hits all 7 sectors}.
    Sector of value v in [0,1) is floor(7 v). 0 in E always hits sector 0.
    Breakpoints: where any e*x crosses a multiple of 1/7, i.e. x = a/(7e)."""
    E = sorted(set(e for e in E if e != 0))
    if not E:
        return Fraction(0)  # only sector 0 ever hit
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        for a in range(0, 7*e+1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(bps)
    tot = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        x = (lo+hi)/2
        sectors = {0}  # speed 0 hits sector 0
        for e in E:
            v = e*x
            v = v - (v.numerator//v.denominator)  # frac part
            sectors.add((v.numerator*7)//v.denominator)
            if len(sectors) == 7: break
        if len(sectors) == 7:
            tot += hi-lo
    return tot

def dyadic_decomp(E):
    """odd-part b -> sorted list of 2-exponents a present in E."""
    d = {}
    for e in E:
        b, a = e, 0
        while b % 2 == 0: b //= 2; a += 1
        d.setdefault(b, []).append(a)
    for b in d: d[b].sort()
    return d

def richness_choose2(E):
    """Sum over odd b of C(chain_len, 2): pairs sharing a doubling chain."""
    d = dyadic_decomp(E)
    s = 0
    for b, exps in d.items():
        L = len(exps); s += L*(L-1)//2
    return s

def richness_maxchain(E):
    """Length of the longest CONTIGUOUS doubling chain b,2b,4b,... in E
    (consecutive exponents). The 'depth' of dyadic structure."""
    d = dyadic_decomp(E)
    best = 0
    for b, exps in d.items():
        # longest run of consecutive integers in exps
        run = 1; mx = 1
        for i in range(1, len(exps)):
            if exps[i] == exps[i-1]+1: run += 1; mx = max(mx, run)
            else: run = 1
        best = max(best, mx)
    return best

def richness_total_chainlen(E):
    """Sum of (chain_len) over odd b with chain_len>=2 minus #chains:
    total 'doubling edges' = number of pairs (e,2e) both in E."""
    Eset = set(E)
    return sum(1 for e in E if 2*e in Eset)

def report(label, E):
    E = sorted(set(E))
    p0 = meas_S7(E)
    print(f"{label:28s} E={E}")
    print(f"   p0=meas(S7)={float(p0):.5f}={p0}")
    print(f"   dyadic={dyadic_decomp(E)}")
    print(f"   richness: C2-pairs={richness_choose2(E)}  doubling-edges(e,2e)={richness_total_chainlen(E)}  maxchain={richness_maxchain(E)}")
    return p0

print("="*78)
print("PART 1: Baselines — consec vs structured alternatives at k=8,9")
print("="*78)
report("consec {1..8}", range(1,9))
report("consec {1..9}", range(1,10))

print()
print("="*78)
print("PART 2: k=8 — exhaustive over ALL 8-subsets of {1..13}: p0 vs richness")
print("="*78)
results = []
universe = list(range(1,14))
for E in itertools.combinations(universe, 8):
    p0 = meas_S7(E)
    results.append((p0, richness_choose2(E), richness_total_chainlen(E), richness_maxchain(E), E))
results.sort(key=lambda t: -t[0])
consec8 = tuple(range(1,9))
print(f"Total 8-subsets of {{1..13}}: {len(results)}")
print(f"\nTOP 6 by p0:")
for p0,r2,re,mc,E in results[:6]:
    tag = " <-- CONSEC" if E==consec8 else ""
    print(f"   p0={float(p0):.5f}  C2={r2:2d} dbl-edges={re} maxchain={mc}  E={E}{tag}")
print(f"\nWhere does CONSEC {consec8} rank?")
for rank,(p0,r2,re,mc,E) in enumerate(results):
    if E==consec8:
        print(f"   CONSEC rank = {rank+1} / {len(results)}  p0={float(p0):.5f}")
        break

# correlation of p0 with each richness measure (Spearman-ish: rank corr sign)
import statistics
def pearson(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy) if sx*sy else 0
p0s=[float(t[0]) for t in results]
print(f"\nPearson corr(p0, richness):")
print(f"   C2-pairs:      {pearson(p0s,[t[1] for t in results]):+.4f}")
print(f"   doubling-edges:{pearson(p0s,[t[2] for t in results]):+.4f}")
print(f"   maxchain:      {pearson(p0s,[t[3] for t in results]):+.4f}")

print()
print("="*78)
print("PART 3: Is the p0-MAX subset the richness-MAX subset? (k=8)")
print("="*78)
maxp0 = results[0][0]
top_p0 = [t for t in results if t[0]==maxp0]
print(f"p0-max value = {float(maxp0):.5f}, achieved by {len(top_p0)} subset(s):")
for p0,r2,re,mc,E in top_p0:
    print(f"   E={E}  C2={r2} dbl-edges={re} maxchain={mc}")
maxr2 = max(t[1] for t in results)
top_r2 = [t for t in results if t[1]==maxr2]
print(f"\nrichness(C2)-max value = {maxr2}, achieved by {len(top_r2)} subset(s); their p0 range:")
r2p0 = sorted(float(t[0]) for t in top_r2)
print(f"   p0 in [{r2p0[0]:.5f}, {r2p0[-1]:.5f}]; is consec among them? {consec8 in [t[4] for t in top_r2]}")
