"""
lrc14_lemmaA_7arc_secondmoment_opus_S185.py  (opus-2026-07-09-S185)

TAIL BOUND for Lemma A via the 7-ARC SECOND MOMENT. Partition the circle into 7 fixed arcs
A_j=[j/7,(j+1)/7). If A_j is empty of {frac(e_i x)} then maxgap>1/7 (bracketing points span a gap >1/7,
STRICT). So with N(x)=#{empty arcs}: nu(E) >= P(N>=1) >= E[N]^2/E[N^2] (Paley-Zygmund). E[N],E[N^2] are
arc-AVOIDANCE moments (covering integrals = theta-sums over the resonance lattice), tractable unlike
maxgap. This computes the PZ bound vs actual nu across the spectrum: does PZ >= nuConsec=0.4425 on the
TAIL (spread>B)? And the decorrelated prediction E[N]~7*(6/7)^13, E[N^2]~7*(6/7)^13+42*(5/7)^13.
"""
import numpy as np
from fractions import Fraction as F
from math import comb

NG=300007; X=(np.arange(NG)+0.5)/NG
NUCONSEC=477/1078

def moments_and_nu(E):
    E=np.array(sorted(set(E)),dtype=float)
    P=(np.outer(E,X))%1.0                 # k x NG  in [0,1)
    arc=np.floor(7*P).astype(int)         # which of 7 arcs each point is in
    # occupancy: for each x, bitmask of occupied arcs
    occ=np.zeros((7,NG),dtype=bool)
    for j in range(7):
        occ[j]= (arc==j).any(axis=0)
    N=7-occ.sum(axis=0)                    # number of EMPTY arcs
    EN=N.mean(); EN2=(N.astype(float)**2).mean()
    PZ = EN*EN/EN2 if EN2>0 else 0.0
    # actual nu
    Ps=np.sort(P,axis=0); mg=np.maximum(np.diff(Ps,axis=0).max(axis=0),(Ps[0]+1)-Ps[-1])
    nu=float((mg>1/7).mean())
    PN1=float((N>=1).mean())
    return EN,EN2,PZ,PN1,nu

# decorrelated prediction
d1=(F(6,7))**13; d2=(F(5,7))**13
EN_dec=7*float(d1); EN2_dec=7*float(d1)+42*float(d2)
print(f"decorrelated: E[N]~7*(6/7)^13={EN_dec:.4f}, E[N^2]~{EN2_dec:.4f}, PZ~{EN_dec**2/EN2_dec:.4f}  (need>={NUCONSEC:.4f})")
print("="*104)
print(f"  {'cluster':>26} {'spread':>6} {'E[N]':>6} {'E[N^2]':>7} {'PZ bound':>8} {'P(N>=1)':>8} {'nu':>7} {'PZ>=nuC?':>8}")
def gap(dims,steps):
    import itertools
    return sorted(set(1+sum(x*s for x,s in zip(c,steps)) for c in itertools.product(*[range(d) for d in dims])))
import random; rng=random.Random(5)
fams=[("consecutive {0..12}",list(range(13))),
      ("near-2-AP",[0,2,3,4,5,6,7,8,9,10,11,12,14]),
      ("{0..11}+20",list(range(12))+[20]),
      ("GAP 4x4",gap([4,4],[1,5])[:13]),
      ("2*AP {0,2..24}",[2*i for i in range(13)]),]
for sp in [20,24,30,40,60,100]:
    fams.append((f"minhunt spread{sp}", None, sp))
for item in fams:
    if len(item)==2:
        name,E=item
    else:
        name,_,sp=item
        # adversarial: minimize PZ at this spread (the worst tail case for our bound)
        best=None;bestPZ=9
        for t in range(40):
            E=sorted(rng.sample(range(0,sp+1),13))
            if max(E)!=sp: continue
            EN,EN2,PZ,PN1,nu=moments_and_nu(E)
            if PZ<bestPZ: bestPZ=PZ;best=(E,EN,EN2,PZ,PN1,nu)
        E,EN,EN2,PZ,PN1,nu=best; 
        print(f"  {name:>26} {max(E):>6} {EN:>6.3f} {EN2:>7.3f} {PZ:>8.4f} {PN1:>8.4f} {nu:>7.4f} {'YES' if PZ>=NUCONSEC else 'no':>8}")
        continue
    EN,EN2,PZ,PN1,nu=moments_and_nu(E)
    print(f"  {name:>26} {max(E):>6} {EN:>6.3f} {EN2:>7.3f} {PZ:>8.4f} {PN1:>8.4f} {nu:>7.4f} {'YES' if PZ>=NUCONSEC else 'no':>8}")
print("="*104)
print("READING: if PZ bound >= nuConsec=0.4425 for spread>B (adversarial minhunt), the 7-arc second moment")
print("PROVES the tail (nu>=PZ>=nuConsec); the crossover B scopes the finite core. Near consecutive PZ<nuC")
print("(expected: PZ<=P(N>=1)<=nu, and nu=nuConsec is the min there) -- that band is the finite core.")
