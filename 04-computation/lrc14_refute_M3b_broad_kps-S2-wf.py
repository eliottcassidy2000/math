"""
Broad DISTINCT-SPEED search for forbidden classes, honestly distinguishing:
  (A) STRUCTURAL witnesses: distinct nonzero residues mod 14 => no tie-break, the
      tournament is fully determined by the loneliness/section geometry. These are the
      strongest refutations.
  (B) TIE-BREAK witnesses: some vote ties resolved by speed order. Still a legitimate
      output of the EXACT map on distinct speeds, but partly an artifact of tie-break.

We search distinct primitive 5-speed sets with speeds in {1..R} for growing R, reference
= every member, and record which forbidden classes appear and of what type. Speeds need
NOT be in {1..13}: the task explicitly allows broader LRC-constrained inputs (more speeds,
larger/smaller, all units). A 5-speed primitive set is a valid LRC instance.
"""
from itertools import combinations
from math import gcd
import random

UNITS=[1,3,5,9,11,13]
def d0(v,a):
    s=(v*a)%14
    return min(s,14-s)
def build(S, ref):
    V=sorted(S); n=len(V)
    w={a:(1 if (ref*a)%14%2==0 else -1) for a in UNITS}
    adj=[[0]*n for _ in range(n)]
    tie_used=False
    for i in range(n):
        for j in range(n):
            if i==j: continue
            x,y=V[i],V[j]; vote=0
            for a in UNITS:
                dx=d0(x,a); dy=d0(y,a)
                if dx>dy: vote+=w[a]
                elif dx<dy: vote-=w[a]
            if vote>0: adj[i][j]=1
            elif vote<0: adj[i][j]=0
            else:
                adj[i][j]=1 if x>y else 0
                tie_used=True
    return adj,V,tie_used
def scoreseq(adj):
    return tuple(sorted(sum(adj[i]) for i in range(len(adj))))
def c3(adj):
    n=len(adj); c=0
    for i in range(n):
        for j in range(n):
            if i==j or not adj[i][j]: continue
            for k in range(n):
                if k in (i,j): continue
                if adj[j][k] and adj[k][i]: c+=1
    return c//3
def H(adj):
    n=len(adj); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            cur=dp[mask][last]
            if not cur: continue
            for nxt in range(n):
                if mask&(1<<nxt): continue
                if adj[last][nxt]: dp[mask|(1<<nxt)][nxt]+=cur
    return sum(dp[full][v] for v in range(n))
def sk(adj): return (H(adj),c3(adj),scoreseq(adj))
def isprim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

FORBIDDEN={(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))}

def classify(S,ref):
    adj,V,tie=build(list(S),ref)
    residues=[v%14 for v in V]
    distinct_nonzero = len(set(residues))==5 and all(r!=0 for r in residues)
    return sk(adj),distinct_nonzero,tie,tuple(residues)

def main():
    structural={f:[] for f in FORBIDDEN}  # distinct nonzero residues, no tie used
    tiebreak={f:[] for f in FORBIDDEN}
    # Exhaustive distinct speeds in {1..R}
    for R in (13,20,27,40):
        cntS=0
        for S in combinations(range(1,R+1),5):
            if not isprim(S): continue
            cntS+=1
            for ref in S:
                key,distinct_nonzero,tie,residues=classify(S,ref)
                if key in FORBIDDEN:
                    rec=(tuple(sorted(S)),ref,residues)
                    if distinct_nonzero and not tie:
                        structural[key].append(rec)
                    else:
                        tiebreak[key].append(rec)
        print(f"R={R}: primitive 5-subsets={cntS}")
        for f in sorted(FORBIDDEN):
            ns=len(structural[f]); nt=len(tiebreak[f])
            print(f"   {f}: structural(distinct-residue,no-tie)={ns}  tie-break/other={nt}")
        print(flush=True)

    print("=== STRUCTURAL witnesses (the strong refutations: distinct nonzero residues, "
          "tournament determined purely by section geometry, NO tie-break) ===")
    for f in sorted(FORBIDDEN):
        if structural[f]:
            print(f"\nFORBIDDEN {f} -> STRUCTURALLY REALIZED. examples:")
            for rec in structural[f][:8]:
                print("   speeds",rec[0],"ref",rec[1],"residues mod14",rec[2])
        else:
            print(f"\nFORBIDDEN {f} -> no structural witness up to R=40 (distinct residues).")

    print("\n=== Forbidden classes realized ONLY via tie-break (weaker) ===")
    for f in sorted(FORBIDDEN):
        if not structural[f] and tiebreak[f]:
            print(f"   {f}: {len(tiebreak[f])} tie-break witnesses, e.g. {tiebreak[f][0]}")
        if not structural[f] and not tiebreak[f]:
            print(f"   {f}: NOT realized at all up to R=40.")

if __name__=="__main__":
    main()
