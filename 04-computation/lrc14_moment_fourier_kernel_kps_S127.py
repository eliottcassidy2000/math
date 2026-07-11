# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont26: WORK THE FOURIER IDENTITY m2 ~ E2 -- honestly.
# (1) Is m2 = E2 an IDENTITY or a proxy?  Derive the exact pair term; test the linear fit R^2.
# (2) RECONCILE my cont.25 corr(F2,E2)=-0.98 (uncontrolled battery) vs opus-S222 corr(L,E2)=-0.58
#     (controlled diam<=13): re-run E2 vs longest-AP head-to-head on CONTROLLED 0-anchored cores.
# (3) The exact pair-collision kernel: does it depend on difference only (=> E2/translation-inv shadow)
#     or on sum too (=> 0-anchored, longest-AP sharp)?  This decides the sharp invariant.
import itertools, math

def sector(y): return min(6, int((y % 1.0) * 7))
def moments(E, N=24001):
    s1 = s2 = 0
    for i in range(N):
        x = (i + 0.5) / N
        n = 6 - len(set(sector(e * x) for e in E) & {1,2,3,4,5,6})
        s1 += n; s2 += n*(n-1)
    m1 = s1/N; m2 = s2/N
    return m1, m2, 4*m1 - m2, 6*m1 - m2          # F2 (kps), L (opus/mac-mini THM-705)
def E2(F):
    d = {}
    for a in F:
        for b in F: d[a+b] = d.get(a+b,0)+1
    return sum(v*v for v in d.values())
def longestAP(F):
    S = set(F); best = 1
    Fl = sorted(F)
    for i in range(len(Fl)):
        for j in range(i+1, len(Fl)):
            d = Fl[j]-Fl[i]; c = 2; nxt = Fl[j]+d
            while nxt in S: c += 1; nxt += d
            best = max(best, c)
    return best
def collide(a, b, N=24001):  # exact inner same-sector collision meas{sector(ax)=sector(bx) in {1..6}}
    c = 0
    for i in range(N):
        x = (i+0.5)/N; sa = sector(a*x); sb = sector(b*x)
        if sa == sb and 1 <= sa <= 6: c += 1
    return c/N
def pear(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/(sx*sy) if sx*sy>0 else 0

def main():
    # (3) pair-collision kernel: is collide(a,b) a function of (b-a) only, or also (a+b)?
    print("(3) PAIR-COLLISION KERNEL collide(a,b) -- difference vs sum dependence:")
    print("     fix diff d=b-a=3, vary a: collide should be CONSTANT if difference-only")
    for a in [1,2,4,7,11,18]:
        print(f"       collide({a},{a+3}) = {collide(a,a+3):.4f}")
    print("     fix a+b=20, vary diff: ")
    for a,b in [(9,11),(8,12),(6,14),(3,17)]:
        print(f"       collide({a},{b}) d={b-a:2d} = {collide(a,b):.4f}")

    # (1)+(2) controlled cores: k=9, 0-anchored within {0..13}
    print("\n(1)/(2) CONTROLLED k=9 cores (0 + 8 from {1..13}), sample:")
    pool = list(itertools.combinations(range(1,14), 8))
    import random; random.seed(1); sample = random.sample(pool, 160)
    cores = [[0]+list(s) for s in sample] + [list(range(9))]
    F2s=[]; Ls=[]; e2s=[]; aps=[]; sss=[]
    # precompute collide table for needed pairs
    coll_cache={}
    for E in cores:
        m1,m2,f2,L = moments(E)
        F2s.append(f2); Ls.append(L); e2s.append(E2(E)); aps.append(longestAP(E))
        ss=0.0
        for a,b in itertools.combinations(E,2):
            key=(a,b)
            if key not in coll_cache: coll_cache[key]=collide(a,b)
            ss+=coll_cache[key]
        sss.append(ss)
    print(f"   n={len(cores)} cores.  Predictors of the residue (consec = argmin F2/L):")
    print(f"     corr(F2, E2)        = {pear(F2s,e2s):+.3f}")
    print(f"     corr(F2, longestAP) = {pear(F2s,[float(a) for a in aps]):+.3f}")
    print(f"     corr(F2, collideSum)= {pear(F2s,sss):+.3f}   <- the EXACT pair term")
    print(f"     corr(L , E2)        = {pear(Ls,e2s):+.3f}   (opus-S222 got -0.58)")
    print(f"     corr(L , longestAP) = {pear(Ls,[float(a) for a in aps]):+.3f}")
    print(f"     corr(L , collideSum)= {pear(Ls,sss):+.3f}")
    # linear fit R^2 of m2-proxy: F2 vs collideSum vs E2
    def r2(xs,ys): return pear(xs,ys)**2
    print(f"\n   R^2(F2 ~ E2) = {r2(F2s,e2s):.3f} ; R^2(F2 ~ collideSum) = {r2(F2s,sss):.3f} ; R^2(F2 ~ longestAP) = {r2(F2s,[float(a) for a in aps]):.3f}")
    print("   VERDICT: the sharp pair invariant is the one with highest R^2 (collideSum = exact Fejer pair term).")

if __name__ == "__main__":
    main()
