from fractions import Fraction as F
from math import gcd
from functools import reduce
def distZ(num,den):
    r=num%den; return F(min(r,den-r),den)
def reach_M(V):
    V=[abs(v) for v in V if v!=0]; n=len(V); dens=set()
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0); bt=None
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=min(distZ(v*m,d) for v in V)
            if mn>best: best=mn; bt=F(m,d)
    return best,bt
def saturated(V): return all(any(v%q==0 for v in V) for q in range(2,15))

thr=F(1,14)
print("=== NEAR-AP SATURATED families: the real LRC(14) hard core? (opus-S131) ===\n")
tests = {
  "{1..13} (AP, NON-sat)": list(range(1,14)),
  "{2..14} consecutive": list(range(2,15)),
  "{3..15} consecutive": list(range(3,16)),
  "{7..19} consecutive": list(range(7,20)),
  "{14..26} consecutive": list(range(14,27)),
  "{1..13} swap 1->14": [14]+list(range(2,14)),
  "{1..12,14}": list(range(1,13))+[14],
  "{1..13,14}\{7}": [i for i in range(1,15) if i!=7],
}
minM=F(1); worst=None
for name,V in tests.items():
    if len(set(V))!=13: print(f"  [skip {name}: {len(set(V))} distinct]"); continue
    g=reduce(gcd,V); Vr=[v//g for v in V]
    M,bt=reach_M(V); sat=saturated(V)
    if sat and M<minM: minM=M; worst=(name,V)
    print(f"  {name:<26} M={str(M):>7}={float(M):.5f}  t={bt}  saturated={sat}  {'<-- M<1/14!' if M<thr else ''}")

# systematic: ALL 13-consecutive {a..a+12} for a=1..40, and AP-with-one-mult14-swap
print("\n  systematic 13-consecutive {a..a+12}:")
for a in range(1,30):
    V=list(range(a,a+13)); M,bt=reach_M(V); sat=saturated(V)
    tag = "SAT" if sat else "non-sat(sieve)"
    print(f"    a={a:>2}: {{{a}..{a+12}}} M={float(M):.5f}={M} t={bt} [{tag}]" + (" <-- BELOW 1/14" if M<thr else ""))
print(f"\n  min M over saturated near-AP tested: {minM}={float(minM):.5f} (thr 1/14={1/14:.5f}); at {worst[0] if worst else '-'}")
