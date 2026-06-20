import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def p0(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; h=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); h.add((v.numerator*7)//v.denominator)
        if len(h)==7: tot+=hi-lo
    return tot
def dyadic_pairs(E):
    """#pairs (e,2e) both in E (the Glaisher doubling content)."""
    S=set(E); return sum(1 for e in S if 2*e in S)
def residues_mod7(E):
    return len({e%7 for e in E})  # how many of the 7 residue classes covered
def prim(E): return reduce(gcd,tuple(E))==1
print("=== k=8: does p_0 correlate with dyadic doubling content / mod-7 residue coverage? ===")
# exhaustive bounded k=8, group by (dyadic pairs, residues mod 7)
rows=[]
for tail in itertools.combinations(range(1,14),7):
    E=(0,)+tail
    if not prim(E): continue
    rows.append((p0(E), dyadic_pairs(E), residues_mod7(E), E))
rows.sort(reverse=True)
consec=tuple(range(8))
print(f"  consec_8: p_0={float(p0(consec)):.5f} dyadicpairs={dyadic_pairs(consec)} resid_mod7={residues_mod7(consec)}")
print("  TOP 6 by p_0:")
for p,d,r,E in rows[:6]:
    print(f"    p_0={float(p):.5f} dyadic={d} resid7={r}  {E}")
print("  does max-p_0 == max-dyadic? max dyadic =",max(r[1] for r in rows),"; consec dyadic=",dyadic_pairs(consec))
# correlation: average p_0 by dyadic count
from collections import defaultdict
byd=defaultdict(list)
for p,d,r,E in rows: byd[d].append(float(p))
print("  avg p_0 by dyadic-pair count:", {d:round(sum(v)/len(v),4) for d,v in sorted(byd.items())})
byr=defaultdict(list)
for p,d,r,E in rows: byr[r].append(float(p))
print("  avg p_0 by #residues-mod-7 covered:", {r:round(sum(v)/len(v),4) for r,v in sorted(byr.items())})
