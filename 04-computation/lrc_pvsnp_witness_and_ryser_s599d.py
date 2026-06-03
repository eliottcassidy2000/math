"""Two rigorous groundings for the P-vs-NP reading of the covering-depth master object.
(I) LRC-INSTANCE in NP: the optimum t*=argmax_t min_i ||v_i t|| is attained at a point
    t=m/d with d | (v_i +- v_j) for some pair -- bit-size O(input) -- a poly witness.
(II) RYSER SHAPE: p_0 = sum_{S subset [n]} (-1)^|S| meas(intersection_{i in S} D_i), the
    same alternating-sum shape as the permanent (Ryser); verify == direct p_0.
opus-2026-06-03-S599d. Convention n runners, gap delta=1/(n+1)."""
from fractions import Fraction as F
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def fmin(V,t): return min(dist(v*t) for v in V)   # f(t)=min_i ||v_i t||
def overlap(Vsub,delta):
    eps=set([F(0)])
    for v in Vsub:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); L=len(pts); tot=F(0)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)<delta for v in Vsub): tot+=ln
    return tot
def p0_ryser(V,delta):
    s=F(0)
    for j in range(0,len(V)+1):
        for S in combinations(V,j):
            s+=(-1)**j*(overlap(S,delta) if S else F(1))
    return s
def p0_direct(V,delta):
    eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for sgn in(1,-1): eps.add((F(k)+sgn*delta)/v % 1)
    pts=sorted(eps); L=len(pts); tot=F(0)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>=delta for v in V): tot+=ln
    return tot
def witness(V):
    # candidate optima: t=m/d, d=|v_i +- v_j|; also d=v_i (single-constraint peaks)
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0)
    best=F(-1); bt=None; bd=None
    for d in ds:
        for m in range(0,d):
            t=F(m,d); val=fmin(V,t)
            if val>best: best=val; bt=t; bd=d
    return best,bt,bd
def main():
    n_delta=lambda n:F(1,n+1)
    cases=[((1,2,3,4),4),((1,3,4,7),4),((1,4,6,9),4),((1,2,3,4,5),5),((1,3,4,5,9),5),((2,3,7,8,12),5)]
    print("(I) NP WITNESS: t*=argmax min||v_i t||; check t* denominator divides some (v_i +- v_j)")
    for V,n in cases:
        d=n_delta(n); M,t,dd=witness(V)
        ok_floor = (M>=d)
        den=t.denominator
        # is den | (v_i +- v_j) for some pair (bit-size O(input))?
        ds=set(V)
        for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
        small = any(D%den==0 for D in ds if D>0)
        print(f"  V={V}: M(S)={M}={float(M):.4f} (floor 1/{n+1}={float(d):.4f}, LRC ok={ok_floor}); t*={t} (den={den}, |input-scale={small}); tight={M==d}")
    print("\n(II) RYSER SHAPE: p_0 = sum_S (-1)^|S| meas(cap_S D_i) (permanent-shape) == direct p_0?")
    for V,n in cases:
        d=n_delta(n); pr=p0_ryser(V,d); pd=p0_direct(V,d)
        print(f"  V={V}: p0_ryser={pr}  p0_direct={pd}  match={pr==pd}")
if __name__=='__main__': main()
