"""
S74b: Is the LRC odd obstruction (the kappa_3 associator) LITERALLY the octonion/Fano associator?
Cayley-Dickson two-towers synthesis: the odd/hard half = loss of ASSOCIATIVITY = the octonion (apex-7).
The Fano plane PG(2,2) has a CYCLIC (Singer) model on Z/7: points = Z/7, the 7 lines = translates of the
QR planar difference set {1,2,4} mod 7. The octonion associator [e_i,e_j,e_k] VANISHES on Fano lines
(collinear = associative) and is +-1 OFF lines (non-collinear = non-associative).
TEST: compute the 3-way joint-emptiness cumulant kappa_3(X_i,X_j,X_k) over the 7 sectors (X_s=1[sector s empty]),
classify triples by Fano-line membership (Singer lines mod 7), and see if the associator structure shows up:
do Fano-LINE triples have systematically different kappa_3 than NON-line triples?
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

def sector_of(p): return int((p%1)*7)
# Singer/Fano lines mod 7 = translates of the QR difference set {1,2,4}
FANO=[frozenset(((1+c)%7,(2+c)%7,(4+c)%7)) for c in range(7)]
def is_fano(tri): return frozenset(tri) in FANO

def emptiness_moments(E):
    """E[prod X_s] for all subsets up to size 3, X_s=1[sector s empty] over the 7 sectors."""
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b)
    p1=np.zeros(7); p2=np.zeros((7,7)); p3={}
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        w=float(x1-x0)
        hit=set(sector_of(e*((x0+x1)/2)) for e in E)
        empty=[s for s in range(7) if s not in hit]
        for s in empty: p1[s]+=w
        for s in empty:
            for u in empty: p2[s][u]+=w
        for tri in combinations(sorted(empty),3): p3[tri]=p3.get(tri,0.0)+w
    return p1,p2,p3

def kappa3(tri,p1,p2,p3):
    i,j,k=tri
    Eijk=p3.get(tuple(sorted(tri)),0.0)
    return (Eijk - p2[i][j]*p1[k] - p2[i][k]*p1[j] - p2[j][k]*p1[i] + 2*p1[i]*p1[j]*p1[k])

print("="*88)
print(" OCTONION/FANO TEST: is kappa_3 (3-way emptiness associator) structured by the Singer/Fano lines?")
print(" Fano (Singer) lines mod 7 = translates of QR {1,2,4}:")
for L in sorted(set(tuple(sorted(l)) for l in FANO)): print("   ",L)
print("="*88)
for name,Eset in [("consec {0..6}",tuple(range(7))),("consec {0..7}",tuple(range(8))),
                  ("consec {0..12}",tuple(range(13)))]:
    p1,p2,p3=emptiness_moments(Eset)
    line_k=[]; nonline_k=[]
    for tri in combinations(range(7),3):
        if 0 in tri: continue            # sector 0 = observer, X_0==0, trivial
        k3=kappa3(tri,p1,p2,p3)
        (line_k if is_fano(tri) else nonline_k).append((tri,k3))
    lm=np.mean([k for _,k in line_k]); nm=np.mean([k for _,k in nonline_k])
    lstd=np.std([k for _,k in line_k]); nstd=np.std([k for _,k in nonline_k])
    print(f"\n {name}:  Fano-line triples (not thru 0): {len(line_k)};  non-line: {len(nonline_k)}")
    print(f"   mean kappa_3  LINE = {lm:+.5f} (sd {lstd:.5f})   NON-LINE = {nm:+.5f} (sd {nstd:.5f})   ratio {lm/nm if nm else float('nan'):+.3f}")
    # the 4 inner Fano lines explicitly
    for tri,k in sorted(line_k): print(f"     LINE {tri}: kappa_3={k:+.5f}")
print("="*88)
print(" If LINE vs NON-LINE kappa_3 differ systematically => the associator carries Fano/octonion structure.")
print(" (apex-7 = octonion imaginary units; 14=2*7 = one Cayley-Dickson doubling past the octonion.)")
