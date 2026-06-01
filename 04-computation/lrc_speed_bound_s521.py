#!/usr/bin/env python3
"""
lrc_speed_bound_s521.py   claudebox-2026-06-01-S521

Investigating the speed bound B(n) -- the external ingredient the multiplicative-
walk methodology needs (see lrc-methodology-needs-a-speed-bound-s521.md).

THREE FINDINGS:

(1) DENOMINATOR-SPEED COUPLING (verified).  The maximal denominator needed for a
    lonely time grows only like log(max speed): the extremal family is
    (1, lcm(2..Q)), speed ~ e^Q, lonely at q ~ Q+1.  So with a speed bound B, the
    methodology's modulus is Q ~ ln B (small).  The bottleneck is B(n) itself, and
    the finite check is over speeds < lcm(2..Q) ~ B.

(2) METHODOLOGY-CHECK = RESIDUES MOD lcm(Q) (proved).  Loneliness at q depends only
    on speeds mod q, hence at all q in Q only on speeds mod L=lcm(Q).  So
    "every primitive system is lonely at some q in Q" is a FINITE check over speed
    sets with entries in [0,L).  Number of sets ~ L^{n-1} ~ B^{n-1}.

(3) THE FAST-RUNNER SPEED BOUND (the real engine).  If the other n-2 movers have a
    common-safe time window I (all simultaneously in [1/n,1-1/n]) and the fastest
    mover v satisfies v*|I| >= 1, then v sweeps the whole circle inside I and is
    safe somewhere in I -> lonely.  So a counterexample needs
        v_max < 1 / |largest common-safe window of the other movers|.
    Empirically |I| ~ 1/(2 v_2nd), giving v_max ~ 2 v_2nd WHEN the window exists.
    GAP: for n>=4 the n-2 slower movers may have EMPTY common window (their bad
    sets can tile the circle), and then this gives no bound -- which is exactly why
    B(n) needs careful casework and is large (super-exponential; n<=7 proven this
    way).  This is the genuine hard core; for n=14 the bound is far beyond
    computation.
"""
from fractions import Fraction as F
from math import gcd, lcm
from functools import reduce

N=4
def safe(v,t,n):
    x=(F(v)*t)%1
    return min(x,1-x)>=F(1,n)
def common_window(speeds,n):
    walls=set([F(0),F(1)])
    for v in speeds:
        for k in range(v+1):
            for s in (F(1,n*v),-F(1,n*v)):
                t=F(k,v)+s
                if 0<=t<=1: walls.add(t)
    W=sorted(walls); best=F(0)
    for a,b in zip(W,W[1:]):
        mid=(a+b)/2
        if all(safe(v,mid,n) for v in speeds): best=max(best,b-a)
    return best

def main():
    print("Speed-bound investigation (claudebox-S521)\n")
    print("(3) Fast-runner bound for n=4: v3 < 1/|I(v1,v2)|, |I|=largest common-safe window.")
    print(f"   {'(v1,v2)':>10} {'|I|':>8} {'v3 must be <':>13}")
    for (v1,v2) in [(1,2),(1,3),(3,5),(7,9),(11,13)]:
        I=common_window([v1,v2],4)
        print(f"   {str((v1,v2)):>10} {str(I):>8} {str(float(1/I) if I>0 else 'unbounded'):>13}")
    print("\n   min largest-window over coprime slow pairs with speeds<=S (=> bound on v3):")
    for S in [5,10,20,40]:
        worst=F(1)
        for v1 in range(1,S+1):
            for v2 in range(v1+1,S+1):
                if gcd(v1,v2)!=1: continue
                I=common_window([v1,v2],4)
                if 0<I<worst: worst=I
        print(f"     slow<= {S:>3}: min|I|={float(worst):.4f}  => v3 < {float(1/worst):.1f}  (~2*S: ratio bound)")
    print("\n   GAP: some slow pairs have EMPTY common window (bad sets tile) -> no bound;")
    print("   that casework is why B(n) is large and the recursion does not close trivially.")
    print("\n(1)+(2): with bound B, modulus Q~ln B, check over speeds<lcm(2..Q)~B (~B^(n-1) sets).")
    print("   So B(n) is the whole bottleneck.  Status: B(n) explicit but super-exponential;")
    print("   LRC proven for n<=7 by such inductive interval arguments; n=14 far beyond reach.")

if __name__ == "__main__":
    main()
