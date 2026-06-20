import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
from importlib import util as _u
_spec=_u.spec_from_file_location("base","04-computation/lrc14_dyadic_richness_test_kps.py")
# reuse functions by re-defining (avoid running its main): just redefine here.

def meas_S7(E):
    E = sorted(set(e for e in E if e != 0))
    if not E: return Fraction(0)
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        for a in range(0, 7*e+1): bps.add(Fraction(a, 7*e))
    bps = sorted(bps); tot = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi==lo: continue
        x=(lo+hi)/2; sectors={0}
        for e in E:
            v=e*x; v=v-(v.numerator//v.denominator); sectors.add((v.numerator*7)//v.denominator)
            if len(sectors)==7: break
        if len(sectors)==7: tot+=hi-lo
    return tot

def dyadic_decomp(E):
    d={}
    for e in E:
        b,a=e,0
        while b%2==0: b//=2; a+=1
        d.setdefault(b,[]).append(a)
    for b in d: d[b].sort()
    return d

print("="*78)
print("WHAT do the dyadic-richness-MAX 8-subsets look like? (these had LOW p0)")
print("="*78)
def c2(E):
    return sum(len(v)*(len(v)-1)//2 for v in dyadic_decomp(E).values())
best=[]; mx=0
for E in itertools.combinations(range(1,14),8):
    r=c2(E)
    if r>mx: mx=r; best=[E]
    elif r==mx: best.append(E)
print(f"C2-max = {mx}; example subsets:")
for E in best[:6]:
    print(f"   E={E} dyadic={dyadic_decomp(E)} p0={float(meas_S7(E)):.5f}")
print("=> These pile MANY speeds onto FEW odd-parts (e.g. b=1: 1,2,4,8; b=3: 3,6,12).")
print("   Few distinct odd-parts => few distinct angular directions => LOW coverage.")

print()
print("="*78)
print("REFINED HYPOTHESIS: p0 favors DISTINCT odd-parts (angular spread) AND")
print("contiguity. Test: p0 vs (#distinct odd parts) and (consecutive-integer run).")
print("="*78)
def n_oddparts(E): return len(dyadic_decomp(E))
def max_consec_run(E):
    E=sorted(set(E)); run=1; mx=1
    for i in range(1,len(E)):
        if E[i]==E[i-1]+1: run+=1; mx=max(mx,run)
        else: run=1
    return mx
def is_consec(E):
    E=sorted(set(E)); return all(E[i]==E[i-1]+1 for i in range(1,len(E)))

rows=[]
for E in itertools.combinations(range(1,14),8):
    rows.append((meas_S7(E), n_oddparts(E), max_consec_run(E), is_consec(E), E))
def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy) if sx*sy else 0
p0s=[float(r[0]) for r in rows]
print(f"Pearson corr(p0, #distinct-odd-parts) = {pearson(p0s,[r[1] for r in rows]):+.4f}")
print(f"Pearson corr(p0, max-consec-run)       = {pearson(p0s,[r[2] for r in rows]):+.4f}")
print(f"Pearson corr(p0, is_consec[0/1])       = {pearson(p0s,[1.0 if r[3] else 0.0 for r in rows]):+.4f}")

print()
print("="*78)
print("DOUBLING-MONOTONICITY: start from consec {1..8}, replace ONE element by an")
print("unrelated value (breaking dyadic chains). Does p0 DECREASE every time?")
print("="*78)
base=list(range(1,9)); base_p0=meas_S7(base)
print(f"base consec {base}: p0={float(base_p0):.5f}")
all_down=True; cnt_down=0; cnt_total=0
for i in range(8):
    for new in range(9,14):
        if new in base: continue
        E=base[:i]+[new]+base[i+1:]
        p=meas_S7(E); cnt_total+=1
        if p<base_p0: cnt_down+=1
        else: all_down=False
print(f"single-element replacements of consec{{1..8}} -> larger value:")
print(f"   {cnt_down}/{cnt_total} STRICTLY decreased p0. All decreased? {all_down}")

print()
print("="*78)
print("DOUBLED-APEX {1..k-1, 2k}: the '2nd-loneliest'. Is it the largest-p0")
print("NON-consec set? Compare to other single-apex perturbations at k=8 (base 1..7).")
print("="*78)
# base {1..7} is k=7 consec; add an 8th element. doubled-apex would be 2*8=16 (out of range)
# Reinterpret: consec is {1..8}; doubled-apex = {1..7, 16}? 16 not allowed (max 13).
# Use the in-range analog: {1..7} u {2*k} with k=8 -> but limited to 13.
# Instead: among all {1..7, x} for x in 8..13, rank by p0:
print("sets {1,..,7, x}, x in 8..13:")
for x in range(8,14):
    E=list(range(1,8))+[x]
    print(f"   x={x:2d}  p0={float(meas_S7(E)):.5f}  dyadic={dyadic_decomp(E)}")
print("(x=8 is consec; x=14 would be doubled-apex 2*7 but exceeds 13.)")
