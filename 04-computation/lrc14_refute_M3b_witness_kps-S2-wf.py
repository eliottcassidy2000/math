"""
Confirm the REFUTATION with a CLEAN witness: distinct primitive speeds, verified
through the EXACT original M3b build_tournament, producing forbidden class
(9,3,(1,1,2,3,3)). Also scan for the most LRC-natural witnesses (small speeds,
primitive, distinct residues so NO tie-break ambiguity).
"""
from itertools import combinations
from math import gcd

UNITS=[1,3,5,9,11,13]
def d0(v,a):
    s=(v*a)%14
    return min(s,14-s)
def build(S, ref):
    V=sorted(S); n=len(V)
    w={a:(1 if (ref*a)%14%2==0 else -1) for a in UNITS}
    adj=[[0]*n for _ in range(n)]
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
            else: adj[i][j]=1 if x>y else 0
    for i in range(n):
        for j in range(i+1,n):
            assert adj[i][j]+adj[j][i]==1
    return adj,V
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

TARGET=(9,3,(1,1,2,3,3))

def main():
    print("Searching DISTINCT-residue primitive 5-speed sets in {1..13} (claim's domain) for TARGET",TARGET)
    found_distinct_res=[]
    found_any=[]
    for S in combinations(range(1,14),5):
        if not isprim(S): continue
        for ref in set(S):
            adj,V=build(S,ref)
            if sk(adj)==TARGET:
                residues=[v%14 for v in V]
                distinct_res = len(set(residues))==5 and all(r!=0 for r in residues)
                rec=(tuple(V),ref,tuple(residues),distinct_res)
                found_any.append(rec)
                if distinct_res: found_distinct_res.append(rec)
    print(f"\nWitnesses within {{1..13}} primitive n=5 (claim's EXACT domain): {len(found_any)}")
    for rec in found_any[:12]:
        print("   speeds",rec[0],"ref",rec[1],"residues",rec[2],"distinct&nonzero:",rec[3])
    print(f"\n  ...of which distinct-nonzero-residue (no tie-break, fully canonical): {len(found_distinct_res)}")
    for rec in found_distinct_res[:12]:
        print("   CLEAN:",rec[0],"ref",rec[1],"residues",rec[2])
    # show one full adjacency + verify H,c3,score
    if found_any:
        S,ref,_,_=found_any[0]
        adj,V=build(list(S),ref)
        print(f"\nFull witness: speeds={V}, ref={ref}")
        print("  adjacency (i->j):")
        for row in adj: print("   ",row)
        print("  H=",H(adj)," c3=",c3(adj)," score=",scoreseq(adj))
        print("  matches forbidden class (9,3,(1,1,2,3,3))? ", sk(adj)==TARGET)

if __name__=="__main__":
    main()
