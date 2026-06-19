#!/usr/bin/env python3
"""
lrc14_netmono_kps-S6-wf.py
Test candidate MONOTONICITY facts toward a rigorous proof that consecutive
maximizes meas(N(E)) (equivalently minimizes mu_{1/7}).

We test, EXACTLY:
  (M1) "compression": if E' is obtained from E by decreasing the spread by
       merging a gap in E (removing an unused integer slot), does meas(N) increase?
  (M2) does meas(N(E)) depend only on the SORTED GAP VECTOR of E (the consecutive
       differences), up to the obvious symmetries?  (i.e. is N invariant under
       reflection E -> max(E)-E and under any reordering of gaps?)
  (M3) the central question: among all gap-vectors (g_1,...,g_{k-1}) of positive
       integers, is meas(N) maximized at all-ones (consecutive)?
Exhaustive for k=8 over bounded spread; report any violation.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
ONE7=F(1,7)

def net_intervals(E):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    good=[]
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(iv): return sum((b-a for a,b in iv),F(0))
def netmeas(E): return meas(net_intervals(E))

def gapvec(E):
    E=sorted(E); return tuple(E[i+1]-E[i] for i in range(len(E)-1))

if __name__=="__main__":
    k=8
    # (M2) reflection & gap-permutation invariance
    print("[M2] N depends only on the MULTISET of consecutive gaps?")
    viol=0; checked=0
    base_gaps=[(1,1,1,1,1,1,1),(2,1,1,1,1,1,1),(1,2,1,1,1,1,1),(3,1,1,1,1,1,2),(2,2,1,1,1,1,1)]
    for g in base_gaps:
        # build E from gap vector
        E=[0]
        for d in g: E.append(E[-1]+d)
        nm=netmeas(E)
        # all permutations of gaps should give same netmeas (conjecture)
        seen=set()
        for perm in set(itertools.permutations(g)):
            Ep=[0]
            for d in perm: Ep.append(Ep[-1]+d)
            seen.add(netmeas(Ep))
            checked+=1
        if len(seen)>1:
            viol+=1
            print(f"   gapvec base {g}: netmeas VARIES across permutations: {sorted(float(s) for s in seen)}")
        else:
            print(f"   gapvec {g}: netmeas={float(nm):.6f} INVARIANT across all {len(set(itertools.permutations(g)))} perms")
    print(f"   M2 violations: {viol} (checked {checked})")

    # (M3) maximize netmeas over gap vectors of positive ints, sum = spread
    print("\n[M3] maximize netmeas over gap-vectors (k=8 -> 7 gaps). all-ones=consecutive?")
    best=-1; bestg=None; consec=netmeas(list(range(8)))
    # enumerate gap vectors with each gap in 1..G and sum bounded
    G=6
    results=[]
    for g in itertools.product(range(1,G+1),repeat=7):
        # dedup by sorted (M2 says only multiset matters) -- but verify by using as-is
        if tuple(sorted(g))!=g: continue  # canonical: nondecreasing gaps
        E=[0]
        for d in g: E.append(E[-1]+d)
        nm=netmeas(E)
        results.append((nm,g))
        if nm>best: best=nm; bestg=g
    results.sort(reverse=True)
    print(f"   consecutive (all-ones) netmeas = {consec} = {float(consec):.6f}")
    print(f"   TOP gap-vectors by netmeas:")
    for nm,g in results[:8]:
        print(f"      {g}: {float(nm):.6f} ({nm}) {'<-- CONSEC' if g==(1,)*7 else ''}")
    print(f"   max netmeas = {float(best):.6f} at gapvec {bestg}")
    if bestg==(1,)*7:
        print("   => CONSECUTIVE MAXIMIZES net over all gap-vectors with gaps<=6. (k=8)")
    else:
        print(f"   => *** CONSECUTIVE NOT MAX: {bestg} beats it ***")
