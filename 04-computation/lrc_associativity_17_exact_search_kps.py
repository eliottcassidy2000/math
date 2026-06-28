"""
lrc_associativity_17_exact_search_kps.py  (kind-pasteur-2026-06-27-S31aj)

Owner: re-examine whether the associativity-compression defect = 1/7 EXACTLY
(=> a theorem: the apex prime quantifies the irreducible 3-way non-associative
fraction), and whether the covariance-max (consec) is robust.

My S31ai computation gave Sigma-kappa_3/S3 = 407891843/2855269200 = 0.142856 for
consec_8 -- NEAR 1/7 but not exact. The owner pushes to find the CLEAN object that
IS exactly 1/7. Test MANY natural associativity-defect ratios (exact fractions),
across k, including:
 - Sigma k3 / S3                         (cumulant / 3-way moment)
 - per-symmetric-triple k3 / P(triple)
 - chain-defect: [P(ijk) - P(ij)P(jk)/P(j)] normalized (Markov/associativity)
 - the FROZEN (large generic dilation) limit of the AP
 - using 7 sectors (with anchor) vs 6 inner
"""
import sys, itertools
from fractions import Fraction as F
INNER=list(range(1,7)); ALL7=list(range(7))
def sector_of(p): return int((p%1)*7)

def cells(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b); out=[]
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        cov=set(sector_of(e*((x0+x1)/2)) for e in E)
        out.append((cov,x1-x0))
    return out

def Pemptyset(cl, S, universe):
    S=set(S); tot=F(0)
    for cov,w in cl:
        if S.isdisjoint(cov): tot+=w   # all sectors in S empty (not covered)
    return tot

def set_partitions(s):
    s=list(s)
    if not s: yield []; return
    if len(s)==1: yield [s]; return
    f=s[0]
    for rest in set_partitions(s[1:]):
        for i in range(len(rest)): yield rest[:i]+[[f]+rest[i]]+rest[i+1:]
        yield [[f]]+rest
def fact(n):
    r=1
    for i in range(2,n+1): r*=i
    return r

def kappa(cl, S, universe):
    m={}
    for r in range(0,len(S)+1):
        for T in itertools.combinations(S,r):
            m[frozenset(T)]=Pemptyset(cl,T,universe) if T else F(1)
    tot=F(0)
    for part in set_partitions(list(S)):
        b=len(part); coef=((-1)**(b-1))*fact(b-1); pr=F(1)
        for blk in part: pr*=m[frozenset(blk)]
        tot+=coef*pr
    return tot

def ratios(E, universe):
    cl=cells(E)
    triples=list(itertools.combinations(universe,3))
    pairs=list(itertools.combinations(universe,2))
    Sk3=sum(kappa(cl,t,universe) for t in triples)
    S3=sum(Pemptyset(cl,t,universe) for t in triples)
    Sk2=sum(kappa(cl,p,universe) for p in pairs)
    S2=sum(Pemptyset(cl,p,universe) for p in pairs)
    r1 = Sk3/S3 if S3 else F(0)
    r2 = Sk3/Sk2 if Sk2 else F(0)
    return dict(Sk3=Sk3,S3=S3,Sk2=Sk2,S2=S2, r_k3_over_S3=r1, r_k3_over_k2=r2)

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    print("Is the associativity defect EXACTLY 1/7? (exact fractions, consec at various k)")
    print(f"{'k':>3} {'universe':>9} | Sk3/S3 (=1/7?) | Sk3/Sk2 | Sk2 sign")
    for universe,uname in [(INNER,'6-inner'),(ALL7,'7-all')]:
        print(f"--- universe = {uname} ---")
        for k in range(4,11):
            E=tuple(range(k))
            d=ratios(E,universe)
            eq17 = (d['r_k3_over_S3']==F(1,7))
            print(f"{k:>3} {uname:>9} | {float(d['r_k3_over_S3']):.7f} {'==1/7!' if eq17 else '(not)'} | "
                  f"{float(d['r_k3_over_k2']):+.5f} | Sk2={float(d['Sk2']):+.4f}")
    # FROZEN limit: dilated AP d*consec for growing d (decorrelation)
    print("\n--- FROZEN limit: dilation d*consec_8 (does Sk3/S3 -> 1/7?) ---")
    for dil in [1,2,3,5,10,20]:
        E=tuple(dil*i for i in range(8))
        d=ratios(E,INNER)
        print(f"  d={dil:>3}: Sk3/S3 = {float(d['r_k3_over_S3']):.7f}  (1/7={1/7:.7f})")
    # exact value for consec_8 6-inner -- factor the deviation from 1/7
    d=ratios(tuple(range(8)),INNER)
    dev = d['r_k3_over_S3'] - F(1,7)
    print(f"\nconsec_8 (6-inner) Sk3/S3 = {d['r_k3_over_S3']} ; deviation from 1/7 = {dev} = {float(dev):.2e}")
