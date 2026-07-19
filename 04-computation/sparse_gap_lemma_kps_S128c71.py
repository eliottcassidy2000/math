#!/usr/bin/env python3
"""sparse_gap_lemma_kps_S128c71.py -- kind-pasteur S128 cont.71.
THE SIX-LINEAR-FUNCTIONS STATEMENT, reduced to a SPARSE-GAP LEMMA.

THE REDUCTION.  Work inside a core-safe component I.  The smallest killer k1 cuts I into
gaps of length EXACTLY 6/(7k1) (its own teeth have width 1/(7k1)).  Suppose some k1-gap
contains at most m teeth of k2,k3,k4.  Those teeth have width <= 1/(7k2) each, so the gap
splits into <= m+1 pieces of total length >= 6/(7k1) - m/(7k2), hence the longest piece is

    L >= ( 6/(7k1) - m/(7k2) ) / (m+1).

The four-comb theorem needs L > 1/(7k4).  With m = 2 this reads

    ( 6/k1 - 2/k2 ) / 3 > 1/k4    <=>    6/k1 - 2/k2 > 3/k4 ,               (*)

which for clustered killers (all ~ k) is 4/k > 3/k -- TRUE with 33% room.  With m = 3 it
reads (6/k1 - 3/k2)/4 > 1/k4, i.e. 3/k > 4/k for clustered -- FALSE.  So EVERYTHING turns
on whether some k1-gap is 2-sparse.

COUNTING HEURISTIC.  Each k1-period holds one tooth of each k_i generically, i.e. 3 foreign
teeth; a fraction ~1/7 of them land inside the k1-TOOTH and are harmless, giving an average
of ~18/7 = 2.571 < 3 per gap.  So a 2-sparse gap should exist.  But for clustered combs the
relative phases are FROZEN (THM-1141), so averaging is not automatic -- it must be checked.
PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def k1_gaps(I,k1):
    """full gaps of comb k1 strictly inside interval I=(A,B)"""
    A,B=I; out=[]
    j=int(A*k1)-1
    while F(j,k1)-F(1,14*k1) < B:
        g0=F(j,k1)+F(1,14*k1); g1=F(j+1,k1)-F(1,14*k1)
        if g0>=A and g1<=B and g1>g0: out.append((g0,g1))
        j+=1
    return out
def foreign_teeth(g,ks):
    """count teeth of the combs in ks that intersect gap g"""
    A,B=g; c=0
    for k in ks:
        j=int(A*k)-1
        while F(j,k)-F(1,14*k) < B:
            x=F(j,k)-F(1,14*k); y=F(j,k)+F(1,14*k)
            if y>A and x<B: c+=1
            j+=1
    return c
random.seed(71)
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
print("### the sparse-gap lemma: does some k1-gap hold <= 2 foreign teeth? ###")
print("  regime            trials  min-over-trials of (min teeth per k1-gap)   #trials with no 2-sparse gap")
for tag,step in [("consecutive",1),("step<=3",3),("step<=8",8),("step<=30",30),("spread x1.3",-1)]:
    mins=[]; bad=0; tot=0
    for _ in range(160):
        P=random.choice(C8); M=max(P); lo=13*M+1
        k1=random.randint(lo,lo+600)
        if step>0:
            ks=[k1]
            for _ in range(3): ks.append(ks[-1]+random.randint(1,step))
        else:
            ks=[k1]
            for _ in range(3): ks.append(int(ks[-1]*13//10)+random.randint(0,5))
        iv=safe_set(P)
        I=max(iv,key=lambda s:s[1]-s[0])
        gs=k1_gaps(I,ks[0])
        if not gs: continue
        tot+=1
        counts=[foreign_teeth(g,ks[1:]) for g in gs]
        mn=min(counts); mins.append(mn)
        if mn>2: bad+=1
    if tot: print("  %-17s %-7d %-43d %d"%(tag,tot,min(mins),bad))
print()
print("### does (*) actually hold numerically?  6/k1 - 2/k2 > 3/k4 ###")
print("  killers                     6/k1-2/k2      3/k4        (*) holds")
for ks in [(157,158,159,160),(371,374,377,379),(550,553,554,558),(157,200,300,600),(157,314,628,1256)]:
    lhs=F(6,ks[0])-F(2,ks[1]); rhs=F(3,ks[3])
    print("  %-27s %-14.8f %-11.8f %s"%(str(ks),float(lhs),float(rhs),lhs>rhs))
print()
print("### the resulting bound vs the truth, on the standing worst case ###")
P=[1,3,5,6,7,8,11,12]; ks=(371,374,377,379)
iv=safe_set(P); I=max(iv,key=lambda s:s[1]-s[0])
gs=k1_gaps(I,ks[0]); counts=[foreign_teeth(g,ks[1:]) for g in gs]
m=min(counts)
bound=(F(6,7*ks[0])-F(m,7*ks[1]))/(m+1)
print("  I = [%s, %s], k1-gaps inside: %d ; teeth per gap: %s"%(I[0],I[1],len(gs),sorted(counts)))
print("  sparsest gap m = %d -> bound L >= %s = %.8f"%(m,bound,float(bound)))
print("  needed 1/(7k4) = %.8f ; 7*k4*bound = %.4f"%(1/(7*ks[3]),float(7*ks[3]*bound)))
print("DONE")
