"""
lrc_k8_associativity_cumulant_tower_kps.py  (kind-pasteur-2026-06-27-S31ai)

The OWNER's "explicit failures of compression of properties beyond commutativity
(associativity)". The empty-sector indicators X_j=1[sector j empty] have a CUMULANT
TOWER whose successive terms are compression failures:
 - kappa_1 (means)        : the 1-body content
 - kappa_2 (covariance)   : COMMUTATIVE pairwise correlation (a+b face; folds to Var, degree-2)
 - kappa_3 (3-way)        : the ASSOCIATOR -- the part of P(i,j,k empty) NOT captured by pairwise;
                            the irreducible 3-cocycle = the ASSOCIATIVITY-compression FAILURE
                            (= the odd -9S3/Worpitzky dip term, the apex-7/octonion residue)
 - kappa_4 (4-way)        : the even +6S4 biquadratic term

Tests for k=8 clusters:
 (1) does consec MAXIMIZE Sigma kappa_2 (total covariance) AND Sigma kappa_3 (total associator)?
 (2) the COMPRESSION DEFECT: how much of S3 is NOT predicted by pairwise (the associativity gap)?
 (3) Newton-Maclaurin: m2=E[(N-3)^2], m4=E[(N-3)^4]; which inequality is TIGHT at consec?
 (4) the structure of kappa_3 over the C(6,3)=20 triples (Fano/octonion pattern?).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb
import numpy as np

INNER=list(range(1,7))  # the 6 inner sectors
def sector_of(p): return int((p%1)*7)

def joint_empty_measures(E):
    """P[S] = meas{x: all sectors in S empty}, for every subset S of INNER (exact)."""
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b)
    P={}
    # accumulate per cell: the set of empty inner sectors
    cells=[]
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        cov=set(sector_of(e*((x0+x1)/2)) for e in E)
        empty=frozenset(s for s in INNER if s not in cov)
        cells.append((empty, x1-x0))
    # P[S] = measure of cells whose empty-set CONTAINS S
    from collections import defaultdict
    # only need up to |S|=4
    def Pof(S):
        S=frozenset(S); tot=F(0)
        for empty,w in cells:
            if S<=empty: tot+=w
        return tot
    return Pof

def cumulants(Pof):
    """joint cumulants of indicators X_j=1[j empty] up to order 4, using P[S]=E[prod_{j in S} X_j]."""
    m={frozenset():F(1)}
    for r in range(1,5):
        for S in itertools.combinations(INNER,r):
            m[frozenset(S)]=Pof(S)
    # joint cumulants via moment-cumulant (set-partition) formula
    def kappa(S):
        S=list(S); n=len(S)
        # kappa = sum over partitions pi of S of (|pi|-1)! (-1)^{|pi|-1} prod_blocks m[block]
        total=F(0)
        for part in set_partitions(S):
            b=len(part); coef=((-1)**(b-1))*factorial(b-1)
            prod=F(1)
            for block in part: prod*=m[frozenset(block)]
            total+=coef*prod
        return total
    k2={}; k3={}; k4={}
    for S in itertools.combinations(INNER,2): k2[S]=kappa(S)
    for S in itertools.combinations(INNER,3): k3[S]=kappa(S)
    for S in itertools.combinations(INNER,4): k4[S]=kappa(S)
    return m,k2,k3,k4

def set_partitions(s):
    s=list(s)
    if len(s)==0: yield []; return
    if len(s)==1: yield [s]; return
    first=s[0]
    for rest in set_partitions(s[1:]):
        for i in range(len(rest)):
            yield rest[:i]+[[first]+rest[i]]+rest[i+1:]
        yield [[first]]+rest
def factorial(n):
    r=1
    for i in range(2,n+1): r*=i
    return r

def Ndist_moments(E):
    Pof=joint_empty_measures(E)
    # N-distribution central moments about 3
    # E[(N-3)^p] via N=sum X_j: easier from q_t
    # reconstruct q from cells
    E2=sorted(set(E)); b=set([F(0),F(1)])
    for e in E2:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        cov=set(sector_of(e*((x0+x1)/2)) for e in E2)
        t=7-len(cov)
        if 0<=t<=6: q[t]+=x1-x0
    qf=[float(x) for x in q]
    m2=sum((t-3)**2*qf[t] for t in range(7))
    m4=sum((t-3)**4*qf[t] for t in range(7))
    return m2,m4,qf

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(31)
    consec=tuple(range(8))
    Pof=joint_empty_measures(consec); m,k2,k3,k4=cumulants(Pof)
    Sk2=float(sum(k2.values())); Sk3=float(sum(k3.values())); Sk4=float(sum(k4.values()))
    m2,m4,qf=Ndist_moments(consec)
    print(f"=== consec_8 cumulant tower (associativity = kappa_3) ===")
    print(f"  Sigma kappa_2 (covariance, COMMUTATIVE) = {Sk2:+.5f}")
    print(f"  Sigma kappa_3 (ASSOCIATOR/3-cocycle)    = {Sk3:+.5f}   <-- associativity-compression failure")
    print(f"  Sigma kappa_4 (4-way)                   = {Sk4:+.5f}")
    print(f"  Newton-Maclaurin: m2=E[(N-3)^2]={m2:.5f}  m4=E[(N-3)^4]={m4:.5f}  m4-m2^2={m4-m2**2:+.5f} (CS slack)")
    # the dip term -9 S3 + 6 S4 in cumulant vs moment: S3=sum_{|S|=3} P[S]
    S3=float(sum(m[frozenset(s)] for s in itertools.combinations(INNER,3)))
    S4=float(sum(m[frozenset(s)] for s in itertools.combinations(INNER,4)))
    print(f"  S3 (3-way joint emptiness, moment)={S3:.5f}  S4={S4:.5f}  dip(-9S3+6S4)={-9*S3+6*S4:+.5f}")
    # associativity DEFECT: S3 vs its pairwise prediction (if associative, k3=0 => S3 predicted by k2,k1)
    print(f"  associativity DEFECT = Sigma kappa_3 / S3 = {Sk3/S3 if S3 else 0:+.4f} (0=associative, !=0 = fails)")
    # does consec maximize Sigma k2 and Sigma k3?
    print(f"\n  Over random k=8 clusters -- consec rank for Sigma-kappa_r:")
    pool=[("consec",consec)]
    for _ in range(250):
        cfg=tuple(sorted([0]+random.sample(range(1,18),7)))
        if len(set(cfg))==8: pool.append(("rand",cfg))
    data=[]
    for nm,E in pool:
        Pf=joint_empty_measures(E); mm,a2,a3,a4=cumulants(Pf)
        data.append((float(sum(a2.values())),float(sum(a3.values())),float(sum(a4.values())),E))
    for idx,lbl in [(0,'Sigma k2 (covariance)'),(1,'Sigma k3 (ASSOCIATOR)'),(2,'Sigma k4')]:
        s=sorted(data,key=lambda d:-d[idx])
        rank=[i for i,d in enumerate(s) if d[3]==consec][0]
        print(f"   {lbl:24s}: consec rank={rank}/{len(data)} (0=max)  consec val={dict((d[3],d[idx]) for d in data)[consec]:+.5f}")
