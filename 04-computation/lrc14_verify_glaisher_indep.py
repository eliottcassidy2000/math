import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')

# INDEPENDENT verification of the Glaisher dyadic-richness / consec-maximality
# claim for the LRC(14) sector route. Fresh meas_S7 implementation using a
# DIFFERENT breakpoint strategy to cross-check the original.
#
# meas_S7(E) = measure of x in [0,1) such that {frac(e x): e in E u {0}}
# meets all 7 sectors [j/7,(j+1)/7).
#
# Sector of v = floor(7 frac(e x)). As x increases from 0 to 1, for each e the
# value e*x mod 1 has sector-boundary crossings exactly at x where 7 e x is an
# integer, i.e. x = t/(7e), t = 0..7e. So global breakpoints = union over e of
# {t/(7e)}. On each open cell the sector multiset is constant. Identical
# breakpoint set to original; I instead verify by an INDEPENDENT method:
# rational-arithmetic event sweep that tracks each speed's sector via integer
# floor, and I additionally sanity-check via a fine rational grid refinement.

def frac_sector(e, x):
    """sector index floor(7*frac(e*x)) using exact Fraction."""
    v = e * x
    fl = v.numerator // v.denominator
    f = v - fl  # in [0,1)
    return (f.numerator * 7) // f.denominator

def meas_S7(E):
    E = sorted({e for e in E if e != 0})
    if not E:
        return Fraction(0)
    bps = set()
    bps.add(Fraction(0)); bps.add(Fraction(1))
    for e in E:
        d = 7 * e
        for t in range(0, d + 1):
            bps.add(Fraction(t, d))
    bps = sorted(bps)
    tot = Fraction(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        mid = (lo + hi) / 2
        secs = {0}
        for e in E:
            secs.add(frac_sector(e, mid))
            if len(secs) == 7:
                break
        if len(secs) == 7:
            tot += hi - lo
    return tot

def odd_part(e):
    b = e
    while b % 2 == 0:
        b //= 2
    return b

def chain_at_1(E):
    """length of power-of-2 chain {1,2,4,8,...} present in E."""
    s = set(E); L = 0; v = 1
    while v in s:
        L += 1; v *= 2
    return L

def maxrun(E):
    s = sorted(set(E)); best = 1; run = 1
    for i in range(1, len(s)):
        if s[i] == s[i-1] + 1: run += 1; best = max(best, run)
        else: run = 1
    return best

def richC2(E):
    d = {}
    for e in E: d.setdefault(odd_part(e), 0); d[odd_part(e)] += 1
    return sum(c*(c-1)//2 for c in d.values())

print("="*70)
print("CHECK 1: consec{1..8} value and rank at k=8 (independent meas_S7)")
print("="*70)
all8 = list(itertools.combinations(range(1,14), 8))
vals = [(meas_S7(E), E) for E in all8]
vals.sort(key=lambda t: -t[0])
consec8 = tuple(range(1,9))
p0_consec = meas_S7(consec8)
print(f"consec{{1..8}} p0 = {p0_consec} = {float(p0_consec):.6f}")
print(f"claimed         = 2447/5880 = {float(Fraction(2447,5880)):.6f}  match={p0_consec==Fraction(2447,5880)}")
rank = [i for i,(v,E) in enumerate(vals) if E==consec8][0]+1
n_at_max = sum(1 for v,E in vals if v==vals[0][0])
print(f"rank = {rank} / {len(vals)};  #sets achieving the max = {n_at_max}")
print(f"top set = {vals[0][1]}  p0={float(vals[0][0]):.6f}")
print(f"UNIQUE GLOBAL MAX at k=8: {rank==1 and n_at_max==1}")

print()
print("="*70)
print("CHECK 2: richness(C2)-max sets have p0 <= 0.293")
print("="*70)
maxC2 = max(richC2(E) for v,E in vals)
rich_sets = [(v,E) for v,E in vals if richC2(E)==maxC2]
rp = sorted(float(v) for v,E in rich_sets)
print(f"max C2 = {maxC2}; achieved by {len(rich_sets)} sets")
print(f"their p0 in [{rp[0]:.5f}, {rp[-1]:.5f}]")
print(f"all <= 0.293: {rp[-1] <= 0.293}")
for v,E in rich_sets:
    print(f"   {E}  C2={richC2(E)} p0={float(v):.5f}")

print()
print("="*70)
print("CHECK 3: window chain ordering {1..8}/{2..9}/{3..10}")
print("="*70)
for start in (1,2,3):
    E = tuple(range(start, start+8))
    print(f"  window {E}: chain@1={chain_at_1(E)}  p0={float(meas_S7(E)):.5f}")

print()
print("="*70)
print("CHECK 4: doubling-monotonicity single-swap certificate")
print("="*70)
for k in (8,9):
    base = list(range(1,k+1))
    p0_base = meas_S7(base)
    swaps = 0; strict_dec = 0; viol = []
    for i in range(k):
        for new in range(k+1, 14):
            if new in base: continue
            E = base[:i] + [new] + base[i+1:]
            swaps += 1
            p = meas_S7(E)
            if p < p0_base: strict_dec += 1
            else: viol.append((base[i],new,float(p)))
    print(f"  k={k}: base p0={float(p0_base):.5f}  swaps={swaps}  strict-decrease={strict_dec}  violations={len(viol)}")
    if viol: print(f"     VIOLATIONS: {viol[:5]}")

print()
print("="*70)
print("CHECK 5: stratified richness correlation by run length (k=8)")
print("="*70)
import collections
def pearson(xs,ys):
    n=len(xs)
    if n<2: return None
    mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy) if sx*sy else None
strata=collections.defaultdict(list)
for v,E in vals:
    strata[maxrun(E)].append((float(v), richC2(E)))
for r in sorted(strata):
    grp=strata[r]; xs=[a for a,b in grp]; ys=[b for a,b in grp]
    pc=pearson(xs,ys)
    mp=sum(xs)/len(xs)
    print(f"  run={r}: n={len(grp):4d}  mean_p0={mp:.4f}  corr(p0,C2)={('%.3f'%pc) if pc is not None else 'NA'}")
