#!/usr/bin/env python3
"""
lrc_n2q_qq_cycle_control_s547.py    oracle-2026-06-01-S547

Vigorously explore HYP-2044's sharp open question: at n=2q (a DOUBLED PRIME), does
the (q,q) cycle type control the loneliest LRC configuration, with apex speed = the
repeated cycle length q?

At n=2q the loneliest config (AP at t*=1/(2q)) is the regular 2q-gon: the 2q points
sit at j/(2q), j=0..2q-1 (the 2q-th roots of unity), observer at 0. The pivot is
q = n/2 = the unique order-2 element of Z_{2q}. We test several precise forms:

 (A) the marked symmetry: the reflection v <-> 2q - v (= position j -> -j). Its TWO
     fixed points should be {0, q} = the OBSERVER and the APEX. So the apex is the
     'co-observer': the unique runner fixed by the same reflection that fixes 0.
 (B) the (q,q) cycle type = rotation-by-2 of the 2q-gon: orbits = even vs odd
     positions, two q-cycles. Is it an automorphism of the (unmarked) half-turn
     relation? And which orbit holds the observer / apex / the bracketing runners?
 (C) the cascade trap: clearing runners 1..2q-1 in order, WHICH runner gets the zero
     conditional clearance (the trapped one)? Is it the apex q?
 (D) the value q: apex speed = n/2 = order-2 element = (q,q) cycle length = Burnside
     gcd(q,q). All equal q because q = n/2 -- 'doubling is pairing', q the pivot.
"""
from fractions import Fraction
from itertools import combinations

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def loneliest_collar(speeds, n):
    thr=Fraction(1,n); W={Fraction(0),Fraction(1)}
    for s in speeds[1:]:
        if s==0: continue
        for m in range(0,abs(s)+1):
            for sg in (1,-1):
                t=Fraction(n*m+sg,n*abs(s))
                if 0<=t<=1: W.add(t)
    Wl=sorted(W); pts=Wl+[(a+b)/2 for a,b in zip(Wl,Wl[1:])]
    best=Fraction(-1); bt=None
    for t in pts:
        c=min((d0(Fraction(s)*t) for s in speeds[1:]),default=Fraction(1))
        if c>best: best=c; bt=t
    return bt,best

def half_turn_adj(positions):
    """positions: dict index->Fraction in [0,1). half-turn i->j iff (p_j-p_i) in (0,1/2); None=tie at 1/2."""
    idx=list(positions); m=len(idx); adj={}
    for a in idx:
        for b in idx:
            if a==b: continue
            f=frac(positions[b]-positions[a])
            adj[(a,b)] = None if f==Fraction(1,2) else (1 if f<Fraction(1,2) else 0)
    return adj

def cascade_trap(speeds, n):
    """clear runners in order; return index (1-based runner) of the first ZERO clearance."""
    runs=speeds[1:]; thr=Fraction(1,n)
    def meas(subset):
        sp=(0,)+tuple(subset); W={Fraction(0),Fraction(1)}
        for s in sp[1:]:
            if s==0: continue
            for m in range(0,abs(s)+1):
                for sg in (1,-1):
                    t=Fraction(n*m+sg,n*abs(s))
                    if 0<=t<=1: W.add(t)
        Wl=sorted(W); tot=Fraction(0)
        for a,b in zip(Wl,Wl[1:]):
            mid=(a+b)/2
            if all(d0(Fraction(s)*mid)>=thr for s in sp[1:]): tot+=(b-a)
        return tot
    prev=Fraction(1); cs=[]
    for i in range(len(runs)):
        cur=meas(runs[:i+1]); c=cur/prev if prev>0 else Fraction(0); cs.append((runs[i],c)); prev=cur
        if c==0: return runs[i], cs
    return None, cs

def main():
    print("n=2q doubled-prime: does the (q,q) cycle type control the loneliest config? (oracle-S547)\n")
    for q in (3,5,7,11):
        n=2*q; speeds=tuple(range(n)); apex=q
        bt,collar=loneliest_collar(speeds,n)
        print(f"=== n={n}=2*{q}: loneliest t*={bt}, collar={collar} (=1/n? {collar==Fraction(1,n)}); apex speed=n/2={apex} ===")
        # positions at t* (the regular 2q-gon)
        pos={j: frac(Fraction(j)*bt) for j in range(n)}
        is_regular = all(pos[j]==Fraction(j,n) for j in range(n))
        print(f"  loneliest config = regular {n}-gon (roots of unity)? {is_regular}")
        # (A) reflection v<->2q-v fixed points
        fixed=[v for v in range(n) if (n-v)%n==v]
        print(f"  (A) reflection v<->{n}-v fixed points = {fixed}  (= observer 0 and APEX {apex}: the apex is the co-observer)")
        # (B) rotation-by-2 = (q,q): orbits
        seen=set(); orbits=[]
        for s in range(n):
            if s in seen: continue
            o=[]; x=s
            while x not in seen: seen.add(x); o.append(x); x=(x+2)%n
            orbits.append(o)
        ctype=tuple(sorted(len(o) for o in orbits))
        evens=[x for x in range(n) if x%2==0]; odds=[x for x in range(n) if x%2==1]
        print(f"  (B) rotation-by-2 cycle type = {ctype} (== (q,q)=({q},{q})? {ctype==(q,q)}); "
              f"observer 0 in {'EVEN' if 0 in evens else 'ODD'} q-cycle, apex {apex} in {'EVEN' if apex in evens else 'ODD'} q-cycle")
        # is rotation-by-2 an automorphism of the unmarked half-turn relation? (check on non-tie pairs)
        adj=half_turn_adj(pos)
        autok=True
        for (a,b),val in adj.items():
            if val is None: continue
            a2,b2=(a+2)%n,(b+2)%n
            if adj.get((a2,b2)) != val: autok=False; break
        print(f"      rotation-by-2 preserves the half-turn relation (an automorphism)? {autok}")
        # nearest bracketing runners of the observer
        nb=[1, n-1]
        print(f"      observer's nearest runners (bracketing) = positions {nb}, parities {['EVEN' if x%2==0 else 'ODD' for x in nb]}")
        # (C) cascade trap
        trap, cs = cascade_trap(speeds, n)
        print(f"  (C) cascade trap (first zero clearance) = runner speed {trap}  (apex={apex}? {trap==apex}); "
              f"is it the LAST runner {n-1}? {trap==n-1}")
        print()
    print("READING below in the reflection.")

if __name__=="__main__":
    main()
