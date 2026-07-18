# opus-2026-07-17-S366 -- HYP-7490: THE STRONGER TOOL -- HIGHER-LEVEL TRUNCATION.
# THM-1060 exhibited primitive families where BONF5 fails (~ -0.17..-0.35) while
# the truth is ~ +0.115.  Full inclusion-exclusion (level 13) is EXACT, hence
# positive, so some odd truncation certifies.  WHICH ONE?
#   uncovered >= BONF_m = 1 - S1 + S2 - ... - S_m   (m odd)
# Computed by DFS over subsets with EMPTY-INTERSECTION PRUNING: if a subset's
# comb intersection is empty, every superset is too.
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
LAM = F(1,14); MODULI=[8,9,10,11,12,13,14]

def teeth01(x):
    w = LAM/x; out=[]
    for j in range(0, x+1):
        a,b = F(j,x)-w, F(j,x)+w
        a,b = max(a,F(0)), min(b,F(1))
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j=[],0,0
    while i<len(u) and j<len(v):
        a,b = max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def meas(iv): return sum(b-a for a,b in iv)

def all_Sk(V):
    """S_k for k=1..n via DFS with empty-intersection pruning."""
    n=len(V); T=[teeth01(v) for v in V]
    S=[F(0)]*(n+1); nodes=[0]
    def rec(i, cur, k):
        for j in range(i, n):
            nxt = inter(cur, T[j]) if k>0 else T[j]
            nodes[0]+=1
            if not nxt: continue            # prune: all supersets empty
            S[k+1] += meas(nxt)
            rec(j+1, nxt, k+1)
    rec(0, None, 0)
    return S, nodes[0]

def uncovered(V):
    live=[(F(0),F(1))]
    for x in V:
        w=LAM/x; out=[]
        for (a,b) in live:
            cur=a
            for j in range(floor((a-w)*x), floor((b+w)*x)+2):
                lo2,hi2=F(j,x)-w, F(j,x)+w
                if hi2<=cur: continue
                if lo2>=b: break
                if lo2>cur: out.append((cur,lo2))
                cur=max(cur,hi2)
                if cur>=b: break
            if cur<b: out.append((cur,b))
        live=out
        if not live: break
    return meas(live)

random.seed(366)
for L in [7, 31, 101, 331]:
    core=[q*L for q in MODULI]; lo,hi=8*L,14*L; extra=[]
    while len(extra)<6:
        x=random.randint(lo,hi)
        if x in core or x in extra or gcd(x,L)!=1: continue
        extra.append(x)
    V=sorted(core+extra)
    S,nodes = all_Sk(V)
    truth = uncovered(V)
    print(f"L={L}: min speed {V[0]}, gcd {reduce(gcd,V)}, "
          f"truth (uncovered) = {float(truth):+.6f}   [{nodes} subset nodes visited]")
    run = F(1); first_pos = None
    for m in range(1, len(V)+1):
        run += (-1)**m * S[m]
        if m % 2 == 1:
            mark = ""
            if run > 0 and first_pos is None:
                first_pos = m; mark = "   <-- FIRST POSITIVE"
            print(f"    BONF{m:2d} = {float(run):+.6f}{mark}")
    print(f"    => the certificate turns positive at level {first_pos}"
          f" (BONF13 = exact = {float(run) if len(V)%2==1 else 0:+.6f})")
    print()
