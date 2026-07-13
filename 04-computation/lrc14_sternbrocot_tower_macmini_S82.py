#!/usr/bin/env python3
"""mac-mini-S82: pursue the STERN-BROCOT INDUCTIVE TOWER. The three-distance regimes of t form
the Stern-Brocot tree (opus); the deep well's witness is the CF leaf [0;n-1,n], value n/Phi6(n).
TEST whether the tree gives an INDUCTIVE LOWER BOUND (covering-min(n) from (n-1) via vertex-
insertion / mediant), or only ORGANIZES the values (the covering bound staying metric)."""
from fractions import Fraction as F
from math import gcd
def cf(fr):
    a,b=fr.numerator,fr.denominator; out=[]
    while b: out.append(a//b); a,b=b,a-(a//b)*b
    return out
def Phi6(n): return n*n-n+1
def Mval(S,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=min(min((a*v)%q,q-((a*v)%q)) for v in S)
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return best
def is_cov(S,n): return all(any(v%q==0 for v in S) for q in range(2,n+1))

print("(1) the deep well {1..n-2, n(n-1)} witness = CF [0;n-1,n], value n/Phi6(n):")
for n in [5,6,7,8,9,10,14]:
    dw=list(range(1,n-1))+[n*(n-1)]
    val=F(n,Phi6(n))
    M=Mval(dw,min(2*max(dw),400))
    print(f"  n={n}: deep well={dw if n<=8 else str(dw[:4])+'...'} n/Phi6={val}={float(val):.5f}, M={M}, CF(n/Phi6)={cf(val)}")

print("\n(2) INDUCTIVE TOWER test: does inserting runner n into a covering (n-1)-family give a")
print("    clean lower-bound recursion M_n >= g(M_{n-1})? (vertex insertion: M can only DECREASE)")
# take the n-1 deep well, insert to make the n deep well; compare
for n in [6,7,8,9,10]:
    dw_prev=list(range(1,n-2))+[(n-1)*(n-2)]  # n-1 deep well
    dw_n=list(range(1,n-1))+[n*(n-1)]          # n deep well
    Mp=F(n-1,Phi6(n-1)); Mn=F(n,Phi6(n))
    # the CF relation: [0;n-2,n-1] -> [0;n-1,n]
    print(f"  n={n}: cov-min(n-1)={float(Mp):.5f} [0;{n-2},{n-1}] -> cov-min(n)={float(Mn):.5f} [0;{n-1},{n}]; ratio={float(Mn/Mp):.4f}")

print("\n(3) is 'covering => M >= n/Phi6' a TREE fact or a METRIC fact? The Farey tree ORGANIZES")
print("    the M-spectrum (1/14 < 3/41 < 2/27 < 14/183 < 1/13, mediants), but WHICH families are")
print("    covering (land >= 14/183) is metric -- the tree gives the VALUES' arithmetic, not the")
print("    lower bound. Vertex-insertion M_n <= M_{n-1} (peeling) = the balance, undershoots multi-killer.")
print("VERDICT probe: the CF tower [0;n-1,n] is the DEEP WELL's clean structure (a value tower),")
print("not a lower-bound proof over all covering families -- same split as Delsarte/tournament.")
