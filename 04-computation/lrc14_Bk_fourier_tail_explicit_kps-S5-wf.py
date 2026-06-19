#!/usr/bin/env python3
"""
LRC(14) B(k): make the FOURIER route's tail EXPLICIT (the named open crux, step 2/3).

From kps-S5 (PROVED, Parseval):
   mu(E) = F(k) + sum_{m in Lambda(E)\0} ghat(m),   Lambda(E)={m in Z^{k-1}: m.e=0}.
F(k)=iid ceiling. The open piece: bound sum_{m!=0} |ghat(m)| over the relation lattice.

This script computes, for the 1-D 'gap indicator' building block, the EXACT Fourier
coefficients and their decay, to test whether an explicit constant C(k) (or a smooth
minorant) gives a CLOSEABLE bound  mu(E) >= F(k) - (tail) > 0  in the large-spread regime,
complementary to the bounded-spread finite check.

We do NOT attempt the full k-dim ghat (hard). Instead we test the SHARPEST usable form:
the smooth-MINORANT route. Build h <= g on T^{k-1} as a product of 1-D minorants of the
'avoid an arc' events won't work (g is a max, not a product). So we test the simplest
rigorous SANDWICH that IS computable: the Beurling-Selberg minorant of a single-arc
avoidance, summed via the averaging identity = Psi. We already showed Psi -> 0 with spread.

CONCLUSION TARGET: confirm that NO single-arc / excess-gap functional is uniform, and
quantify the decay rate, so the verdict 'covering-system angle alone cannot close B(k)'
is rigorous and the residual is precisely the bulk-equidistribution (Fourier-tail) crux.

Also: directly test the 'two-regime' decomposition numerically --
   regime LARGE-SPREAD (relation-free): mu -> F(k); how fast? (Erdos-Turan rate)
   regime BOUNDED-SPREAD: finite check (mac-mini).
The boundary is the SHORT-RELATION patterns of Lambda(E).
"""
from fractions import Fraction as F
from math import gcd
G0=F(2,7)
def _frac(q): return q-q.__floor__()
def _collision_bps(E):
    E=sorted(set(E)); k=len(E); bp=set([F(0),F(1)])
    for i in range(k):
        for j in range(i+1,k):
            d=E[j]-E[i]
            if d==0: continue
            for m in range(0,d+1): bp.add(F(m,d))
    return bp
def _gap_eq_bps(E):
    E=sorted(set(E)); bp=set(); diffs=set()
    for i in range(len(E)):
        for j in range(len(E)):
            if i!=j: diffs.add(abs(E[j]-E[i]))
    for D in diffs:
        if D==0: continue
        n=0
        while True:
            cand=[(F(n)+F(2,7))/D,(F(n)+F(5,7))/D]
            for x in cand:
                if F(0)<=x<F(1): bp.add(x)
            if min(cand)>=F(1): break
            n+=1
            if n>D+2: break
    return bp
def gaps_at(E,x):
    pts=sorted(set(_frac(e*x) for e in E))
    if len(pts)==1: return [F(1)]
    g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1]); return g
def mu_exact(E):
    bp=sorted(_collision_bps(E)|_gap_eq_bps(E)); bp=[b for b in bp if F(0)<=b<F(1)]
    tot=F(0); pts=bp+[F(1)]
    for a,b in zip(bp,pts[1:]):
        if b<=a: continue
        if max(gaps_at(E,(a+b)/2))>G0: tot+=(b-a)
    return tot
def normalize(E):
    E=sorted(set(E)); g=0
    for e in E: g=gcd(g,e)
    return [e//g for e in E] if g else E
def header(t): print("\n"+"="*74); print(t); print("="*74)

# exact F(k)
from math import comb
def iid_floor(k):
    s=F(0)
    for j in range(k+1):
        base=1-j*G0
        if base<=0: break
        s+=(-1)**j*comb(k,j)*base**(k-1)
    return 1-s

if __name__=="__main__":
    header("F(k) exact (the large-spread ceiling/limit)")
    for k in range(3,14):
        print(f"  k={k:2d}: F(k)={iid_floor(k)} = {float(iid_floor(k)):.6f}")

    header("LARGE-SPREAD CONVERGENCE: relation-free E -> mu -> F(k). Rate?")
    # relation-free shapes of increasing spread; measure |mu - F(k)|.
    import random
    random.seed(11)
    for k in (4,5,6):
        Fk=iid_floor(k)
        print(f"  k={k}, F(k)={float(Fk):.6f}:")
        for scale in (10,50,200,1000):
            # build a relation-free-ish set: 0,1, then powers-ish/coprime spread
            E=[0,1]
            v=1
            while len(E)<k:
                v=v*scale+random.randint(1,scale)
                E.append(v)
            E=normalize(E)
            m=mu_exact(E)
            print(f"     scale~{scale:5d}: E(maxE={E[-1]:8d}) mu={float(m):.6f}  |mu-F|={float(abs(m-Fk)):.6f}")

    header("MIN over BOUNDED-SPREAD vs F(k): is the bounded-spread min the true inf?")
    # the known bounded-spread minimizers (perforated near-APs). For each k, scan a
    # modest perforation family and report the min, compare to F(k).
    from itertools import combinations
    for k in range(4,11):
        Fk=iid_floor(k)
        cap=k+5
        best=None
        for sub in combinations(range(1,cap+1),k-1):
            E=normalize([0]+list(sub))
            if len(E)<k: continue
            m=mu_exact(E)
            if best is None or m<best[1]: best=(E,m)
        print(f"  k={k:2d}: bounded-spread min mu={float(best[1]):.6f}={best[1]}  vs F(k)={float(Fk):.6f}"
              f"   min<F? {'YES (bounded-spread is the binding inf)' if best[1]<Fk else 'no'}")
        print(f"          argmin E={best[0]}")
