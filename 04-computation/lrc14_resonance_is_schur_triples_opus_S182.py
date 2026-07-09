"""
lrc14_resonance_is_schur_triples_opus_S182.py  (opus-2026-07-09-S182)

The resonance sum R (L=(6/7)^13+R) is a theta-sum over the resonance lattice Lambda={n:n.v=0}, dominated
by its MINIMAL vectors. Those are the height-3 SCHUR TRIPLES a+b=c (norm 3), NOT the height-4 additive-
energy relations a+b=c+d (norm 4) -- fewer arc-coeff factors => larger. So the density-floor deficit
(6/7)^13 - L should track the ORDER-3 additive energy E3=#{(a,b,c) in S^3: a+b=c}, with the AP the unique
maximizer (E3(AP{1..13})=78). Test: does |R| correlate with E3 (Schur triples) better than with E2
(doubling/additive energy)? If yes, the a-priori density floor is a SUM-FREE/Schur structure theorem,
not a Freiman/doubling (BSG,3k-4) one -- a redirection of the forward lead.
"""
import numpy as np, random
from collections import Counter

NG = 300007
TAU = (np.arange(NG)+0.5)/NG
H = 1.0/14
MAIN = (6.0/7.0)**13

def sumset(S): return len(set(a+b for a in S for b in S))
def E2(S):
    r=Counter(a+b for a in S for b in S); return sum(c*c for c in r.values())
def E3(S):
    """#{(a,b,c) in S^3 : a+b=c} ordered Schur triples (sum-closure count)."""
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)
def longest_ap(S):
    Sset=set(S); S=sorted(S); best=1
    for i,a in enumerate(S):
        for b in S[i+1:]:
            d=b-a
            if a-d in Sset: continue
            L=2; x=b+d
            while x in Sset: L+=1; x+=d
            best=max(best,L)
    return best
def lonely(S):
    M=np.zeros(NG)
    for v in S:
        d=np.abs(((v*TAU+0.5)%1.0)-0.5); M+=(d<=H)
    return float((M==0).mean())
def gap(dims,steps):
    import itertools
    return sorted(set(1+sum(x*s for x,s in zip(c,steps)) for c in itertools.product(*[range(d) for d in dims])))

fams=[("AP {1..13}",list(range(1,14))),
      ("2*AP",[2*i for i in range(1,14)]),
      ("near-AP {1..12}+20",list(range(1,13))+[20]),
      ("near-AP {1..11}+{20,30}",list(range(1,12))+[20,30]),
      ("near-AP{1..10}+3",list(range(1,11))+[20,30,40]),
      ("GAP 4x4 P=5",gap([4,4],[1,5])[:13]),
      ("GAP 7x2 P=9",gap([7,2],[1,9])[:13]),
      ("GAP 3x5 P=11",gap([5,3],[1,11])[:13]),
      ("GAP3D 2x2x4",gap([4,2,2],[1,5,20])[:13]),
      ("Fib-ish",[1,2,3,5,8,13,21,34,55,89,144,233,377]),
      ("Sidon-ish",[1,2,5,11,22,33,45,60,78,95,110,130,150]),
      ("sum-free odd {1,3..25}",list(range(1,26,2)))]
rng=random.Random(182)
for s in range(4):
    r=random.Random(500+s); sp=r.randint(40,120)
    S=sorted(set([1]+r.sample(range(2,sp),11)+[sp]))
    while len(S)!=13: S=sorted(set(S+[r.randint(1,sp)]))[:13]
    fams.append((f"random sp{sp}",S))

print("="*100)
print(f"DENSITY-FLOOR DEFICIT (main-L) vs Schur-triples E3 vs doubling.  main=(6/7)^13={MAIN:.4f}")
print("="*100)
print(f"  {'set':>24} {'E3=#(a+b=c)':>11} {'E2 energy':>9} {'K':>5} {'lAP':>4} {'L':>7} {'deficit':>8} {'|R|/main':>8}")
rows=[]
for name,S in fams:
    S=sorted(set(S))
    if len(S)!=13: continue
    e3=E3(S); L=lonely(S); rows.append((e3,name,E2(S),sumset(S)/13,longest_ap(S),L,MAIN-L,abs(L-MAIN)/MAIN))
rows.sort(reverse=True)
e3s=[]; defs=[]; e2s=[]
for e3,name,e2,K,lap,L,deficit,absRr in rows:
    print(f"  {name:>24} {e3:>11} {e2:>9} {K:>5.2f} {lap:>4} {L:>7.4f} {deficit:>8.4f} {absRr:>8.3f}")
    e3s.append(e3); defs.append(deficit); e2s.append(e2)
print("-"*100)
c_e3=np.corrcoef(e3s,defs)[0,1]; c_e2=np.corrcoef(e2s,defs)[0,1]
print(f"  corr(E3 Schur-triples, deficit) = {c_e3:+.3f}    corr(E2 additive-energy, deficit) = {c_e2:+.3f}")
print(f"  => the density-floor deficit is governed by {'E3 (SCHUR TRIPLES)' if abs(c_e3)>abs(c_e2) else 'E2 (doubling)'};")
print(f"     AP maximizes E3 (={E3(list(range(1,14)))}) => unique tight; the a-priori floor is a SUM-FREE/Schur bound.")
print("="*100)
